import re
from itertools import chain

"""
This file contains functions that parse Siemens protocol headers (in ASCONV
and XPROT format).

Acknowledgement:
Most of the string parsing code is taken from William Clarke's pymapvbvd
project (https://github.com/wexeee/pymapvbvd).
"""


def parse_twix_hdr(file):
    import re
    import struct
    import numpy as np
    _, n_buffer = np.fromfile(file, dtype=np.uint32, count=2)
    xprotocol = dict()
    pattern = br'(\w{4,})\x00(.{4})'
    pos = file.tell()
    for _ in range(n_buffer):
        tmp = file.read(48)
        matches = re.search(pattern, tmp, re.DOTALL)
        name = matches.group(1).decode('latin1')
        buf_len = struct.unpack('<I', matches.group(2))[0]
        pos += len(matches.group())
        file.seek(pos)
        buf = file.read(buf_len).decode('latin1')
        xprotocol[name] = parse_buffer(buf)
        xprotocol["{}_raw".format(name)] = buf
        pos += buf_len

    return xprotocol


def try_cast(value, key):
    if key.startswith('t'):
        try:
            value = value.strip('"')
        except ValueError:
            pass
    elif key.startswith('b'):
        try:
            value = bool(value)
        except ValueError:
            pass
    elif key.startswith('l') or key.startswith('ul'):
        try:
            value = int(value)
        except ValueError:
            pass
    elif key.startswith('uc'):
        try:
            if value.startswith('0x'):
                value = int(value, 16)
            else:
                value = int(value)
        except ValueError:
            pass
    else:  # try to convert everything else to float
        # obsolete: elif key.startswith('d') or key.startswith('fl'):
        try:
            value = float(value)
        except ValueError:
            pass
    return value


def update_ascconv(prot, key, value, last_string=None):
    if '__attribute__' in key:
        return
    if len(key) > 1:
        if isinstance(key[0], int):  # int -> list
            while len(prot) < key[0]+1:
                if isinstance(key[1], int):
                    prot.append(list())
                else:
                    prot.append(dict())
            update_ascconv(prot[key[0]], key[1:], value, last_string)
        else:  # string -> dict
            last_string = key[0]
            if key[0] not in prot:
                if isinstance(key[1], int):
                    prot[key[0]] = list()
                else:
                    prot[key[0]] = dict()
            update_ascconv(prot[key[0]], key[1:], value, last_string)
    else:
        if isinstance(key[0], int):
            while len(prot) < key[0]+1:
                prot.append(list())
        else:
            last_string = key[0]

        # remove a (for array) from string
        if last_string.startswith('a'):
            last_string = last_string[1:]

        prot[key[0]] = try_cast(value, last_string)


def parse_ascconv(buffer):
    vararray = re.finditer(r'(?P<name>\S*)\s*=\s*(?P<value>\S*)\n', buffer)
    mrprot = dict()
    for v in vararray:
        # now split array name and index (if present)
        vvarray = re.finditer(r'(?P<name>\w+)(\[(?P<ix>[0-9]+)\])?',
                              v.group('name'))
        currKey = []
        for vv in vvarray:
            currKey.append(vv.group('name'))
            if vv.group('ix') is not None:
                currKey.append(int(vv.group('ix')))

        if len(currKey) > 0:
            update_ascconv(mrprot, currKey, v.group('value'))

    return mrprot


def parse_xprot(buffer):
    xprot = {}
    tokens = re.finditer(
        r'<Param(?:Bool|Long|String)\."(\w+)">\s*{([^}]*)', buffer)
    tokensDouble = re.finditer(
        r'<ParamDouble\."(\w+)">\s*{\s*(<Precision>\s*[0-9]*)?\s*([^}]*)',
        buffer)
    # tokensArray = re.finditer(
    #     r'<ParamArray\."(\w+)">\s+{\s+<Visible>.\"(true|false)\"\s+<MinSize>.(\d+)\s+<MaxSize>.(\d+)\s+<Default>.<(\w+).\"(\w*)\">',
    #     buffer)
    # for t in tokensArray:
    #     print(t.groups())
    alltokens = chain(tokens, tokensDouble)

    for t in alltokens:
        name = t.group(1)

        value = re.sub(r'("*)|( *<\w*> *[^\n]*)', '', t.groups()[-1])
        value = re.sub(r'[\t\n\r\f\v]*', '',
                       value.strip())

        if name.startswith('a'):
            out = list()
            for v in value.split():
                out.append(try_cast(v, name[1:]))
            value = out
        else:
            value = try_cast(value, name)

        xprot.update({name: value})

    return xprot


def parse_buffer(buffer):
    reASCCONV = re.compile(
        r'### ASCCONV BEGIN[^\n]*\n(.*)\s### ASCCONV END ###', re.DOTALL)

    ascconv = reASCCONV.search(buffer)
    if ascconv is not None:
        prot = parse_ascconv(ascconv.group(0))
    else:
        prot = dict()

    xprot = reASCCONV.split(buffer)
    if xprot is not None:
        xprot = ''.join([found for found in xprot])
        prot2 = parse_xprot(xprot)

        prot.update(prot2)

    return prot
