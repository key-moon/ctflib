from itertools import accumulate, product
import random
from typing import Any, Generator, TypeGuard
from sage.all import matrix, vector, GF, Integer
from typing import overload

def is_int(val: Any) -> TypeGuard[int]:
  return isinstance(val, (int, Integer))

def int_to_vec(val: int, bits: int):
  assert int(val).bit_length() <= bits
  return vector(GF(2), [val >> i & 1 for i in range(bits)])
def vec_to_int(vec: vector):
  return sum([int(item) << i for i, item in enumerate(vec)])

Details = tuple[int, ...] | list[int]

def form_vals(val_or_vals: "list[LinearFunc | int] | LinearFunc | int", shape: Details):
  vals = val_or_vals if isinstance(val_or_vals, list) else [val_or_vals]
  has_affine = shape[0] == 1
  expected_len = len(shape) - (1 if has_affine else 0)
  assert len(vals) == expected_len
  if has_affine: vals.insert(0, 1)
  return vals

def form_vec(int_vals: list[int], bits_list: Details):
  return vector(sum([list(int_to_vec(val, bits)) for val, bits in zip(int_vals, bits_list)], start=[]))

class LinearFunc:
  def __init__(self, m, in_bits_details: Details):
    self.m = m
    self.in_bits_details = in_bits_details
    assert self.in_bits == sum(self.in_bits_details)
  @property
  def has_affine(self): return self.in_bits_details[0] == 1
  def _ensure_affine(self): assert self.has_affine
  def _with_same_details(self, m):
    return LinearFunc(m, self.in_bits_details)
  def const(self, val: int, val_bits: int):
    self._ensure_affine()
    const = int_to_vec(val, val_bits).column()
    return LinearFunc(const.augment(matrix.zero(GF(2), val_bits, self.in_bits - 1)), self.in_bits_details)
  def _cast(self, val: "LinearFunc | int"):
    if isinstance(val, LinearFunc):
      assert self.shape == val.shape
      return val
    return self.const(val, self.out_bits)
  @property
  def in_bits(self): return self.m.ncols()
  @property
  def out_bits(self): return self.m.nrows()
  @property
  def shape(self):
    return (self.in_bits_details, self.out_bits)

  # binary
  def __xor__(self, other: "LinearFunc | int"): # self ^ other
    other = self._cast(other)
    assert self.shape == other.shape
    return self._with_same_details(self.m + other.m)
  def __lshift__(self, width: int): # self << width
    width = min(width, self.out_bits)
    return self._with_same_details(matrix.zero(GF(2), width, self.in_bits).stack(self.m[:-width]))
  def __rshift__(self, width: int): # self >> width
    width = min(width, self.out_bits)
    return self._with_same_details(self.m[width:].stack(matrix.zero(GF(2), width, self.in_bits)))
  def _or_and_helper(self, value: int, replace_on: int):
    replace_row = vector(GF(2), self.in_bits)
    replace_row[0] = replace_on
    return self._with_same_details(matrix([replace_row if bit == replace_on else row for row, bit in zip(self.m, int_to_vec(value, self.out_bits))]))
  def __or__(self, value: int) -> "LinearFunc":
    return self._or_and_helper(value, 1)
  def __and__(self, value: int) -> "LinearFunc":
    return self._or_and_helper(value, 0)
  def __invert__(self):
    return self ^ ((1 << self.out_bits) - 1)
  def __not__(self):
    return ~self
  def rotr(self, width: int):
    width %= self.out_bits
    return self._with_same_details(self.m[width:].stack(self.m[:width]))
  def rotl(self, width: int):
    return self.rotr(-width)
  
  # gf
  def gfmul(self, mul: int, mod: int):
    assert int(mod).bit_length() == self.out_bits + 1
    x = matrix.zero(GF(2), 1, self.out_bits - 1).stack(matrix.identity(GF(2), self.out_bits - 1)).augment(int_to_vec(mod, self.out_bits + 1)[:-1].column())
    cur = matrix.identity(GF(2), self.out_bits)
    result = matrix.zero(GF(2), self.out_bits, self.in_bits)
    for i in range(self.out_bits):
      if mul >> i & 1:
        result += cur * self.m
      cur *= x
    return self._with_same_details(result)

  # bit selection
  def __getitem__(self, index: int | slice | tuple[bool, ...] | tuple[int, ...] | list[bool] | list[int]):
    if is_int(index):
      assert 0 <= index < self.out_bits
      return self._with_same_details(self.m[index].row())
    elif isinstance(index, slice):
      return self._with_same_details(self.m[index])
    elif isinstance(index, (tuple, list)):
      if all(is_int(elem) for elem in index):
        m = matrix([self.m[i] for i in index])
      else:
        assert len(index) == self.out_bits
        m = matrix([row for row, select in zip(self.m, index) if select])
      return self._with_same_details(m)
    else:
      raise TypeError("Index must be an int, slice, tuple, or list")
  def select(self, zero: int, one: int, bits: int):
    assert int(zero).bit_length() <= bits
    assert int(one).bit_length() <= bits
    self._ensure_affine()
    zero_vector = vector(GF(2), [0] * self.in_bits)
    one_vector = vector(GF(2), [1] + [0] * (self.in_bits - 1))
    rows = []
    for i in range(bits):
      zero_bit, one_bit = zero >> i & 1, one >> i & 1
      if zero_bit and one_bit: rows.append(one_vector)
      if zero_bit and not one_bit: rows.append(self.m[i] ^ one_vector)
      if not zero_bit and one_bit: rows.append(self.m[i])
      if not zero_bit and not one_bit: rows.append(zero_vector)
    return matrix(rows)

  # evaluation / composition
  def __call__(self, *val_or_vals: "LinearFunc | int") -> Any:
    vals = form_vals(list(val_or_vals), self.in_bits_details)
    if all([not isinstance(val, LinearFunc) for val in vals]):
      vec = form_vec(vals, self.in_bits_details) # type: ignore
      return vec_to_int(self.m * vec)
    else:
      mat = matrix(GF(2), 0, self.in_bits)
      in_context = next(filter(lambda x: isinstance(x, LinearFunc), vals))
      assert isinstance(in_context, LinearFunc)
      for val, val_bit in zip(vals, self.in_bits_details):
        c = val if isinstance(val, LinearFunc) else in_context.const(val, val_bit)
        assert c.in_bits_details == in_context.in_bits_details
        assert c.out_bits == val_bit
        mat = mat.stack(c.m)
      return LinearFunc(self.m * mat, in_context.in_bits_details)
  def __mul__(self, val: "LinearFunc"):
    return val(self)
  def __pow__(self, times: int):
    if times == 0: return self.const(1, self.out_bits)
    if times == 1: return self
    if times % 2 == 0:
      return (self * self) ** (times // 2)
    return self * (self ** (times - 1))

  def __repr__(self):
    inp = ", ".join("1" if i == 0 and bits == 1 else f"BV({bits})" for i, bits in enumerate(self.in_bits_details))
    oup = f"BV({self.out_bits})"
    return f"({inp}) -> {oup}"

# TODO: 複数解ある場合にエラーを吐く / iterate するオプション(all)
@overload
def interpolate_function(input_or_inputs: list[list[int]], input_bit_or_bits: Details, outputs: list[int], output_bit: int, affine: bool=False): ...
@overload
def interpolate_function(input_or_inputs: list[int], input_bit_or_bits: int, outputs: list[int], output_bit: int, affine: bool=False): ...
def interpolate_function(input_or_inputs, input_bit_or_bits, outputs: list[int], output_bit: int, affine=False):
  if is_int(input_or_inputs[0]) and is_int(input_bit_or_bits):
    inputs: list[list[int]] = [[input] for input in input_or_inputs]
    input_bits: Details = (input_bit_or_bits,)
  else:
    inputs: list[list[int]] = input_or_inputs
    input_bits: Details = input_bit_or_bits
  assert len(input_or_inputs) == len(outputs)
  if affine:
    inputs = [[1, *input] for input in inputs]
    input_bits = (1, *input_bits)
  input_mat = matrix([form_vec(input, input_bits) for input in inputs]).T
  output_mat = matrix([int_to_vec(output, output_bit) for output in outputs]).T
  # find X s.t. X * in = out
  mat = input_mat.solve_left(output_mat)
  return LinearFunc(mat, input_bits)

@overload
def inverse(function_or_functions: list[LinearFunc], output_or_outputs: list[int], all=False) -> list[int]: ...
@overload
def inverse(function_or_functions: LinearFunc, output_or_outputs: int, all=False) -> list[int]: ...
@overload
def inverse(function_or_functions: list[LinearFunc], output_or_outputs: list[int], all=True) -> Generator[list[int], Any, list[int] | None]: ...
@overload
def inverse(function_or_functions: LinearFunc, output_or_outputs: int, all=True) -> Generator[list[int], Any, list[int] | None]: ...

def inverse(function_or_functions, output_or_outputs, all=False) -> list[int] | Generator[list[int], Any, list[int] | None]:
  if isinstance(function_or_functions, LinearFunc) and is_int(output_or_outputs):
    functions: list[LinearFunc] = [function_or_functions]
    outputs: list[int] = [output_or_outputs]
  else:
    functions: list[LinearFunc] = function_or_functions # type: ignore
    outputs: list[int] = output_or_outputs
  shape, in_bits, has_affine = functions[0].shape, functions[0].in_bits, functions[0].has_affine
  if has_affine:
    out_shape = [1]
    mat = matrix(GF(2), [[1] + [0] * (in_bits - 1)])
    out_vec = form_vec([1, *outputs], [1, *[f.out_bits for f in functions]])
  else:
    out_shape = []
    mat = matrix(GF(2), 0, in_bits)
    out_vec = form_vec(outputs, [f.out_bits for f in functions])
  for func in functions:
    assert func.shape == shape
    out_shape.append(func.out_bits)
    mat = mat.stack(func.m)
  inp = mat.solve_right(out_vec)
  inds = list(accumulate(functions[0].in_bits_details))
  if not has_affine: inds.insert(0, 0)

  hom = mat.right_kernel()
  basis = hom.basis()
  if all:
    for coeffs in product(GF(2), repeat=len(basis)):
      combo = sum(c * b for c, b in zip(coeffs, basis))
      new_inp = inp + combo
      yield [vec_to_int(new_inp[start:stop]) for start, stop in zip(inds[:-1], inds[1:])]
  else:
    if len(basis) != 0:
      print("[!] linear equation has non unique solution. put all=True to get all solutions")
    return [vec_to_int(inp[start:stop]) for start, stop in zip(inds[:-1], inds[1:])]

def get_variables(*bits: int, affine=True):
  in_bits_details = (1, *bits) if affine else bits
  variables = []
  prev_bits = 0
  for i, bit in enumerate(in_bits_details):
    if not (affine and i == 0):
      m = matrix.block([[matrix.identity(GF(2), bit) if i == j else matrix.zero(GF(2), bit, block_bit) for j, block_bit in enumerate(in_bits_details)]])
      variables.append(LinearFunc(m, in_bits_details))
    prev_bits += bit
  return tuple(variables)

def get_variable(bits: int, affine=True):
  return get_variables(bits, affine=affine)[0]

if __name__ == "__main__":
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

  ground_truth = xorshift(0xdeadbeefcafebabe)
  assert xorshift(var)(0xdeadbeefcafebabe) == ground_truth
  assert xorshift_2(xorshift_1(0xdeadbeefcafebabe)) == ground_truth
  assert (xorshift_2(xorshift_1))(0xdeadbeefcafebabe) == ground_truth
  assert f(0xdeadbeefcafebabe) == ground_truth

  state = 0xdeadbeefcafebabe
  for _ in range(100):
    state = xorshift(state)
  
  assert ((xorshift_2(xorshift_1)) ** 100)(0xdeadbeefcafebabe) == state
  assert inverse(xorshift_2(xorshift_1) ** 100, state)[0] == 0xdeadbeefcafebabe

  print("OK!")
