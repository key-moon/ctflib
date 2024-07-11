from setuptools import setup, find_packages

setup(
    name="keymoon-ctf-libs",
    version="0.0.1",
    install_requires=["pwntools", "ptrlib", "pycryptodome"],
    extras_require={
    },
    packages=find_packages(),
    entry_points={
    }
)
