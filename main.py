import numpy

from AESTypes import AESTypes
from galois_fields import ExtendedGaloisField
from polinomials import Polynomials

class AES:
    def __init__(self, key_len):
        self.key_len = key_len
        # self.source_data_block = numpy.random.randint(low=2, size=128)
        # self.key = numpy.random.randint(low=2, size=key_len)
        self.source_data_block = self.hex_to_bin('3243f6a8885a308d313198a2e0370734')
        self.key = self.hex_to_bin('2b7e151628aed2a6abf7158809cf4f3c')
        self.round_keys = []
        self.gf = ExtendedGaloisField(prime_order=2, n=pow(2, 8) - 1, primitive_poly=0b100011101)
        self.poly = Polynomials(gf=self.gf)

    def hex_to_bin(self, str):
        bin_list = []
        int_list = [int(i, 16) for i in str]
        for i in range(len(int_list)):
            x = [int(d) for d in bin(int_list[i])[2:].zfill(4)]
            bin_list += x

        return bin_list

    def generate_state_matrix(self):

        matrix = []

        for i in range(AESTypes.NB):
            tmp_list = []
            for j in range(AESTypes.NB*8*i, AESTypes.NB*8*(i+1), 8):
                cell_value = int("".join(str(x) for x in self.source_data_block[j:j+8]), 2)
                tmp_list.append(cell_value)

            matrix.append(tmp_list)

        return matrix


    def sub_bytes(self, dword):
        sbox_res = []
        for i in range(len(dword)):
            sbox_res.append(AESTypes.Sbox[dword[i]])

        return sbox_res

    def inv_sub_bytes(self, dword):
        sbox_res = []
        for i in range(len(dword)):
            sbox_res.append(AESTypes.InvSbox[dword[i]])

        return sbox_res

    def xor_dword(self, dword1, dword2):
        xor_res = []
        for i in range(len(dword1)):
            xor_res.append(dword1[i] ^ dword2[i])

        return xor_res

    def shift_word(self, dword, steps):

        for i in range(steps):
            tmp = dword[0]
            dword.pop(0)
            dword.append(tmp)
        return dword

    def inv_shift_word(self, dword, steps):

        for i in range(steps):
            tmp = dword[-1]
            dword.pop()
            dword.insert(0, tmp)
        return dword

    def shift_rows(self, state_matrix):

        shifted_state_matrix = state_matrix

        for i in range(len(state_matrix)):
            tmp_row = []
            for j in range(len(state_matrix)):
                tmp_row.append(state_matrix[j][i])

            tmp_row = self.shift_word(tmp_row, i)

            for j in range(len(tmp_row)):
                shifted_state_matrix[j][i] = tmp_row[j]

        return shifted_state_matrix

    def inv_shift_rows(self, state_matrix):

        shifted_state_matrix = state_matrix

        for i in range(len(state_matrix)):
            tmp_row = []
            for j in range(len(state_matrix)):
                tmp_row.append(state_matrix[j][i])

            tmp_row = self.inv_shift_word(tmp_row, i)

            for j in range(len(tmp_row)):
                shifted_state_matrix[j][i] = tmp_row[j]

        return shifted_state_matrix

    def mult_in_gf(self, a, b):
        a = [int(n) for n in bin(a)[2:]]
        b = [int(n) for n in bin(b)[2:]]

        mult = self.poly.bin_poly_mult(a, b)
        mult.reverse()

        mod_mult = self.poly.bin_poly_div(mult, [1, 1, 0, 1, 1, 0, 0, 0, 1])
        mod_mult[1].reverse()

        return int("".join(str(x) for x in mod_mult[1]), 2)


    def mix_columns(self, state_matrix):

        mixed_matrix = []

        for i in range(len(state_matrix)):
            tmp_column = []
            for coefs in AESTypes.POLY_COEFs:
                s_i = self.mult_in_gf(coefs[0], state_matrix[i][0]) ^ self.mult_in_gf(coefs[1], state_matrix[i][1]) ^ \
                      self.mult_in_gf(coefs[2], state_matrix[i][2]) ^ self.mult_in_gf(coefs[3], state_matrix[i][3])

                tmp_column.append(s_i)
            mixed_matrix.append(tmp_column)

        return mixed_matrix


    def inv_mix_columns(self, state_matrix):

        mixed_matrix = []

        for i in range(len(state_matrix)):
            tmp_column = []
            for coefs in AESTypes.INV_POLY_COEFs:
                s_i = self.mult_in_gf(coefs[0], state_matrix[i][0]) ^ self.mult_in_gf(coefs[1], state_matrix[i][1]) ^ \
                      self.mult_in_gf(coefs[2], state_matrix[i][2]) ^ self.mult_in_gf(coefs[3], state_matrix[i][3])

                tmp_column.append(s_i)
            mixed_matrix.append(tmp_column)

        return mixed_matrix

    def key_expansion(self):
        nk = AESTypes.NK[self.key_len]

        for i in range(AESTypes.NK[self.key_len]):
            tmp = []
            for j in range(i*32, (i+1)*32, 8):
                tmp.append(int("".join(str(x) for x in self.key[j:j + 8]), 2))
            self.round_keys.append(tmp)

        i = nk
        while i < AESTypes.NB * (AESTypes.KEY_ROUND[self.key_len] + 1):
            tmp = self.round_keys[i - 1].copy()
            if i % nk == 0:
                tmp = self.shift_word(tmp, 1)
                tmp = self.sub_bytes(tmp)
                tmp = self.xor_dword(tmp, AESTypes.Rcon[int(i / nk)])
            elif (nk > 6) and (i % nk == 4):
                tmp = self.sub_bytes(tmp)
            self.round_keys.append(self.xor_dword(self.round_keys[i-nk], tmp))
            i += 1

        return

    def add_round_key(self, state_matrix, round_key):

        round_key_matrix = []

        for i in range(AESTypes.NB):
            row_xor = self.xor_dword(state_matrix[i], round_key[i])
            round_key_matrix.append(row_xor)

        return round_key_matrix

    def run_round(self, rk_matrix):
        matrix = rk_matrix.copy()
        for i in range(AESTypes.KEY_ROUND[self.key_len]):
            tmp_matrix = []
            for j in range(len(matrix)):
                tmp_matrix.append(self.sub_bytes(matrix[j]))
            matrix = tmp_matrix
            print("{} SubBytes:".format(i+1))
            print_matrix(matrix)

            matrix = self.shift_rows(matrix)
            print("{} ShiftRows: ".format(i+1))
            print_matrix(matrix)

            if i < AESTypes.KEY_ROUND[self.key_len] - 1:
                matrix = self.mix_columns(matrix)
                print("{} MixColumns:".format(i+1))
                print_matrix(matrix)

            matrix = self.add_round_key(matrix, self.round_keys[AESTypes.NB * (i+1):AESTypes.NB * (i+2)])
            print("{} Key: ".format(i + 1))
            print_matrix(self.round_keys[AESTypes.NB * (i+1):AESTypes.NB * (i+2)])
            print("{} AddRoundKey: ".format(i+1))
            print_matrix(matrix)

        return matrix

    def inv_round(self, res_matrix):

        matrix = res_matrix.copy()

        for i in range(AESTypes.KEY_ROUND[self.key_len]-1, -1, -1):
            matrix = self.add_round_key(matrix, self.round_keys[AESTypes.NB * (i+1):AESTypes.NB * (i+2)])
            print("{} AddRoundKey: ".format(i + 1))
            print_matrix(matrix)

            if i != AESTypes.KEY_ROUND[self.key_len]-1:
                matrix = self.inv_mix_columns(matrix)
                print("{} MixColumns:".format(i + 1))
                print_matrix(matrix)

            matrix = self.inv_shift_rows(matrix)
            print("{} ShiftRows: ".format(i + 1))
            print_matrix(matrix)

            tmp_matrix = []
            for j in range(len(matrix)):
                tmp_matrix.append(self.inv_sub_bytes(matrix[j]))
            matrix = tmp_matrix
            print("{} SubBytes:".format(i + 1))
            print_matrix(matrix)

        matrix = self.add_round_key(matrix, self.round_keys[0:4])

        return matrix


def print_matrix(matrix):
    for i in range(4):
        print(str(hex(matrix[0][i])) + " " + str(hex(matrix[1][i])) + " " + str(hex(matrix[2][i])) +
              " " + str(hex(matrix[3][i])))

def encrypt(aes):
    aes.key_expansion()
    print("Key 0: ")
    print_matrix(aes.round_keys)

    state_matrix = aes.generate_state_matrix()
    print("State matrix: ")
    print_matrix(state_matrix)

    rk_matrix = aes.add_round_key(state_matrix, aes.round_keys[:AESTypes.NB])
    print("First AddRoundKey: ")
    print_matrix(rk_matrix)
    res_matrix = aes.run_round(rk_matrix=rk_matrix)

    return res_matrix

def decrypt(aes, res_matrix):

    matrix = aes.inv_round(res_matrix)
    print("Result:")
    print_matrix(matrix)
    return

def main():
    key_len = 128
    aes = AES(key_len=key_len)

    # print("Message: 3243f6a8885a308d313198a2e0370734")
    # print("Key: 2b7e151628aed2a6abf7158809cf4f3c")

    res_matrix = encrypt(aes)
    decrypt(aes, res_matrix)



    return


if __name__ == "__main__":
    main()