from setuptools import setup

setup(
    name="gerenuk-reloaded",
    version="0.2.0",
    author="Jeet Sukumaran, Morzh_Tulenev",
    author_email="jeetsukumaran@gmail.com",
    packages=["gerenuk"],
    scripts=[
        "bin/gerenuk-simulate.py",
        "bin/gerenuk-reject.py",
    ],
    url="https://github.com/yourname/gerenuk-reloaded",
    license="MIT",
    description="Fork of gerenuk",
    long_description=open("README.txt", encoding="utf-8").read(),
    long_description_content_type="text/plain",
    python_requires=">=3.9",
    install_requires=[],
    extras_require={
        "test": ["pytest>=7"],
    },
)
