import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--ids", nargs=3)

args = parser.parse_args()

print repr(args)
