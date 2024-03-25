from setuptools import setup, find_packages

# set long description to readme file 
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="msrDynamics",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    description="Python API to JiTCDDE for neutronic and thermal hydraulic systems",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Luke Labrie-Cleary",
    author_email="lplcleary@gmail.com",
    license="MIT",
    include_package_data=True,
    install_requires=[
        'jitcdde',
        'numpy',
        'chspy'
    ],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
        "Environment :: Console",
    ],
    python_requires='>=3.6',
)
