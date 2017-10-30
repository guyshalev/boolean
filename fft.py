import numpy, math
maj = [1,1,1,-1,1,-1,-1,-1]


def calc_fft(table):
    assert 2 ** int(math.log(len(table), 2)) == len(table)
    atable = numpy.array(table)
    atable = numpy.resize(atable, [2] * int(math.log(len(table), 2)))
    return numpy.fft.fftn(atable) / len(table)


def address_func(chooser, options):
    chooser = int('0b' + ''.join(map(lambda x: str(int(0.5* (1-x))), chooser)),2)
    return options[chooser]
    

def majority(inp):
    x = sum(inp)
    return int(x / abs(x))


def tribes(inp, s = 3, w = 3):
    assert len(inp) == s * w
    for i in xrange(s):
        if all([x==-1 for x in inp[i*w :(i+1)*w]]):
            return -1
    return 1


print calc_fft(maj)
