#-*- coding: iso-8859-1 -*-
import sys


def main():
    filename = sys.argv[1]
    parameter = sys.argv[2]
    value = sys.argv[3]

    lines = []

    f = open(filename, 'r')
    parameternotfound = True
    for line in f:
        if len(line.split()) > 1:
            if line.split()[0] == r'$' + parameter:
                parameternotfound = False
                # replace
                if value == 'delete':
                    line = ''
                else:
                    linearr = line.split()
                    linearr[1] = value + "\n"
                    line = linearr[0] + " " + linearr[1]
        lines.append(line)
    if parameternotfound:
        # append parameter
        lines.append(r'$' + parameter + ' ' + value + '\n')
    f.close()
    f = open(filename, 'w')

    for line in lines:
        f.write(line)

    f.close()


if __name__=='__main__':
    main()



