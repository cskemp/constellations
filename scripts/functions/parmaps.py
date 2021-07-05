import numpy as np

pm = {}
pm['bmap'] = lambda bs: [2 if 4*b < 2 else np.min([4*b, 5.0]) for b in bs]
# b2 for single values rather than lists
pm['b2map'] = lambda b: 2 if 4*b < 2 else np.min([4*b, 5.0])
pm['rmap'] = lambda r: (5 ** r -1)
pm['r_tqmap'] = lambda r: (5 ** r -1)
pm['nmap'] = lambda n: int(200 * n + 0.00000001)
pm['dmap'] = lambda d: 30 * d
# CODE parameters
pm['tmap'] = lambda t: t if t < 1 and t > 0.4 else 0.4
pm['bfmap'] = lambda t: t
pm['hmap'] = lambda t: t
# kNN parameters
pm['kmap'] = lambda t: t

invpm = {}
invpm['bmap'] = lambda b: 0.25 * b
invpm['b2map'] = lambda b: 0.25 * b
invpm['rmap'] = lambda r: np.log(r + 1) / np.log(5)
invpm['r_tqmap'] = lambda r: np.log(r + 1) / np.log(5)
invpm['nmap'] = lambda n: n / 200
invpm['dmap'] = lambda d: d / 30
# CODE parameters
invpm['tmap'] = lambda t: t
invpm['bfmap'] = lambda t: t
invpm['hmap'] = lambda t: t
# kNN parameters
invpm['kmap'] = lambda t: t

def parmaps():
    return pm, invpm




