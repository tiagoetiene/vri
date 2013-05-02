#! /usr/bin/env python2.7
from math import *
import plot
import sys
from time import sleep


def main():

    x = []
    y = []
    while True:
        sleep(1)
        l = sys.stdin.readline()
        data = [float(d) for d in l.split()]
        if len(x) == 0:
            x = data
        elif len(y) == 0:
            y = data
        else:
            break
    print x
    print y
    plot.plot(x, y)
    

if __name__ == "__main__":
    main()
