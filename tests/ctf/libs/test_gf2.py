import pytest
import random
from sage.all import matrix, vector, GF, Integer
from ctf.libs.gf2 import (
    int_to_vec, vec_to_int, LinearFunc, interpolate_function,
    inverse, get_variables, get_variable
)

print("Test file is being executed!")

class TestConversionFunctions:
    def test_int_to_vec(self):
        # Test conversion of integers to vectors
        assert list(int_to_vec(0, 4)) == [0, 0, 0, 0]
        assert list(int_to_vec(1, 4)) == [1, 0, 0, 0]
        assert list(int_to_vec(5, 4)) == [1, 0, 1, 0]  # 5 = 0b0101
        assert list(int_to_vec(15, 4)) == [1, 1, 1, 1]  # 15 = 0b1111
        
        # Test with larger numbers
        assert list(int_to_vec(0xA, 8)) == [0, 1, 0, 1, 0, 0, 0, 0]  # 10 = 0b00001010
        
        # Test with bit length exactly matching bits parameter
        assert list(int_to_vec(15, 4)) == [1, 1, 1, 1]
        
        # Test with assertion error when bit length exceeds bits parameter
        with pytest.raises(AssertionError):
            int_to_vec(16, 4)  # 16 = 0b10000, which is 5 bits
    
    def test_vec_to_int(self):
        # Test conversion of vectors to integers
        assert vec_to_int(vector(GF(2), [0, 0, 0, 0])) == 0
        assert vec_to_int(vector(GF(2), [1, 0, 0, 0])) == 1
        assert vec_to_int(vector(GF(2), [1, 0, 1, 0])) == 5
        assert vec_to_int(vector(GF(2), [1, 1, 1, 1])) == 15
        
        # Test with larger vectors
        assert vec_to_int(vector(GF(2), [0, 1, 0, 1, 0, 0, 0, 0])) == 10
    
    def test_roundtrip_conversion(self):
        # Test that converting from int to vec and back gives the original int
        for i in range(100):
            num = random.randint(0, 2**16 - 1)
            bits = num.bit_length()
            assert vec_to_int(int_to_vec(num, bits)) == num


class TestLinearFunc:
    def test_init(self):
        # Test initialization with valid parameters
        m = matrix(GF(2), [[1, 0], [0, 1]])
        func = LinearFunc(m, (2,))
        assert func.in_bits == 2
        assert func.out_bits == 2
        assert func.in_bits_details == (2,)
        
        # Test initialization with affine term
        m = matrix(GF(2), [[1, 0, 0], [0, 1, 0]])
        func = LinearFunc(m, (1, 2))
        assert func.in_bits == 3
        assert func.out_bits == 2
        assert func.in_bits_details == (1, 2)
        assert func.has_affine == True
        
        # Test initialization with assertion error when in_bits doesn't match sum of in_bits_details
        with pytest.raises(AssertionError):
            LinearFunc(m, (1, 1))
    
    def test_const(self):
        # Test const method
        var = get_variable(8, affine=True)
        const = var.const(5, 8)
        assert const(0) == 5
        assert const(1) == 5
        assert const(255) == 5
    
    def test_binary_operations(self):
        # Test XOR operation
        var = get_variable(8, affine=True)
        result = var ^ 5
        assert result(0) == 5
        assert result(1) == 4
        assert result(5) == 0
        
        # Test shift operations
        var = get_variable(8, affine=True)
        left_shift = var << 2
        right_shift = var >> 2
        assert left_shift(5) == 20
        assert right_shift(20) == 5
        
        # Test bitwise operations
        var = get_variable(8, affine=True)
        and_result = var & 0x0F
        or_result = var | 0xF0
        not_result = ~var
        assert and_result(0xFF) == 0x0F
        assert or_result(0x0F) == 0xFF
        assert not_result(0x55) == 0xAA
        
        # Test rotation operations
        var = get_variable(8, affine=True)
        rotr_result = var.rotr(2)
        rotl_result = var.rotl(2)
        assert rotr_result(0b11110000) == 0b00111100
        assert rotl_result(0b00111100) == 0b11110000
    
    def test_composition(self):
        # Test function composition
        var = get_variable(8, affine=True)
        f1 = var ^ 5
        f2 = var << 1
        composed = f2(f1)
        assert composed(10) == ((10 ^ 5) << 1)
        
        # Test function power
        var = get_variable(8, affine=True)
        f = (var << 1) ^ (var >> 1)
        f_squared = f ** 2
        assert f_squared(10) == f(f(10))


class TestInterpolation:
    def test_interpolate_function(self):
        # Define a simple function to interpolate
        def simple_func(x):
            return (x ^ (x << 1)) & 0xFF
        
        # Generate input-output pairs
        inputs = [random.randint(0, 255) for _ in range(20)]
        outputs = [simple_func(x) for x in inputs]
        
        # Interpolate the function
        f = interpolate_function(inputs, 8, outputs, 8)
        
        # Test the interpolated function on the original inputs
        for x, y in zip(inputs, outputs):
            assert f(x) == y
        
        # Test the interpolated function on new inputs
        for _ in range(10):
            x = random.randint(0, 255)
            assert f(x) == simple_func(x)
    
    def test_interpolate_function_with_affine(self):
        # Define a simple affine function to interpolate
        def affine_func(x):
            return (x + 10) & 0xFF
        
        # Generate input-output pairs
        inputs = [random.randint(0, 255) for _ in range(20)]
        outputs = [affine_func(x) for x in inputs]
        
        # Interpolate the function with affine=True
        f = interpolate_function(inputs, 8, outputs, 8, affine=True)
        
        # Test the interpolated function on the original inputs
        for x, y in zip(inputs, outputs):
            assert f(x) == y
        
        # Test the interpolated function on new inputs
        for _ in range(10):
            x = random.randint(0, 255)
            assert f(x) == affine_func(x)


class TestInverse:
    def test_inverse_simple(self):
        # Create a simple linear function
        var = get_variable(8, affine=True)
        f = var ^ 5
        
        # Test inverse
        for x in range(10):
            y = f(x)
            x_inv = inverse(f, y)[0]
            assert x_inv == x
    
    def test_inverse_complex(self):
        # Create a more complex function
        var = get_variable(8, affine=True)
        f = (var ^ (var << 1)) & 0xFF
        
        # Test inverse
        for x in range(10):
            y = f(x)
            x_inv = inverse(f, y)[0]
            assert f(x_inv) == y
    
    def test_inverse_multiple_solutions(self):
        # Create a function with multiple solutions
        var = get_variable(8, affine=True)
        f = var & 0xFE  # Clear the least significant bit
        
        # Test inverse with all=True
        y = f(42)
        solutions = list(inverse(f, y, all=True))
        assert len(solutions) == 2
        assert 42 in [sol[0] for sol in solutions]
        assert 43 in [sol[0] for sol in solutions]
        assert all(f(sol[0]) == y for sol in solutions)


class TestVariables:
    def test_get_variable(self):
        # Test get_variable with affine=True
        var = get_variable(8, affine=True)
        assert var.in_bits == 9
        assert var.out_bits == 8
        assert var.has_affine == True
        
        # Test get_variable with affine=False
        var = get_variable(8, affine=False)
        assert var.in_bits == 8
        assert var.out_bits == 8
        assert var.has_affine == False
        
        # Test that the variable acts as identity function
        for x in range(10):
            assert var(x) == x
    
    def test_get_variables(self):
        # Test get_variables with affine=True
        var1, var2 = get_variables(8, 16, affine=True)
        assert var1.in_bits == 25  # 1 + 8 + 16
        assert var1.out_bits == 8
        assert var2.in_bits == 25
        assert var2.out_bits == 16
        
        # Test that the variables extract the correct parts
        for x in range(5):
            for y in range(5):
                assert var1(x, y) == x
                assert var2(x, y) == y


class TestExamples:
    def test_xorshift_example(self):
        # Test the example from the module
        def xorshift(x):
            x ^= (x << 7) & 0xffff_ffff_ffff_ffff
            x ^= x >> 9
            return x
        
        var = get_variable(64)
        
        xorshift_1 = var ^ (var << 7)
        xorshift_2 = var ^ (var >> 9)
        
        inputs = [random.randrange(0, 2**64) for _ in range(70)]
        outputs = list(map(xorshift, inputs))
        
        f = interpolate_function(inputs, 64, outputs, 64)
        
        test_input = 0xdeadbeefcafebabe
        ground_truth = xorshift(test_input)
        
        assert xorshift(var)(test_input) == ground_truth
        assert xorshift_2(xorshift_1(test_input)) == ground_truth
        assert (xorshift_2(xorshift_1))(test_input) == ground_truth
        assert f(test_input) == ground_truth
        
        # Test multiple iterations
        state = test_input
        for _ in range(10):  # Reduced from 100 to 10 for faster tests
            state = xorshift(state)
        
        assert ((xorshift_2(xorshift_1)) ** 10)(test_input) == state
        assert inverse(xorshift_2(xorshift_1) ** 10, state)[0] == test_input
