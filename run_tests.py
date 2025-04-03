import sys
import os
import pytest

# Add the src directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

# Run the tests
result = pytest.main(["-v", "tests/ctf/libs/test_gf2.py"])
print(f"Test result: {result}")
