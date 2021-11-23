#!/usr/bin/python
# did not use in the end!
import sys
# COMMAND: python FormatBED_selected.py 'BED: selected bed with no overlapping regions' 'Out_BED: header: Chr, Start, End, CN_estimation'
def main ():
    for line in sys.stdin:
        #tmpLine = line
        tmpRead = line.split('\t')
        tmpRead[12] = tmpRead[12].split(':')[0]
        print(tmpRead)
main()

# def main():
#     tmpRead = []
#     tmpLine = []
#     tmpCounter = 0
#     f = open(sys.argv[1])
#     line = f.readline()
#     prevline = line
#     tmpRead = line.split('\t')
#     #test reading files:
#     #line = f.readline()
#     #currContent = line.split('\t')
#     print(tmpRead[1])
#     f.close()
#     return
# main()
