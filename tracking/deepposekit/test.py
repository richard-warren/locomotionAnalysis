import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--video_name')
parser.add_argument('--model_name')
args = parser.parse_args()

print(args.video_name)

##

if 1:
    print('yes')