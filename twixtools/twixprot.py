def twixprot(f, hdr_len, version_is_ve):
    """
    Parses xprotocol of siemens twix-file object and returns a protocol dict.

    Currently only (the first occurence of) the ascconv-protocol
    data is read. Seems to be enough information for now.
    """
    # output dict:
    ascconv = dict()

    # find begin of ascconv data
    search_str = '### ASCCONV BEGIN'
    if version_is_ve:
        search_str += ' object=MrProtDataImpl@MrProtocolData'

    start_pos = f.tell()
    while str(f.readline()).find(search_str) == -1:
        if f.tell()-start_pos > hdr_len:
            return ascconv

    for line in f:
        if f.tell()-start_pos > hdr_len:
            break
        line = line.decode('utf-8').rstrip()
        if line.find('### ASCCONV END') > -1:
            break

        cConv = ascconv
        parts = line.split()
        splitparts = parts[0].split('.')
        for i, var in enumerate(splitparts):
            if var.startswith('a'):
                tmp = var.strip(']').split('[')
                key = tmp[0]
                if len(tmp) > 1:
                    index = int(tmp[1])
                else:
                    index = 0
                if key not in cConv.keys():
                    cConv.update({key: list()})
                if i < len(splitparts) - 1:
                    for k in range(len(cConv[key]) - 1, index):
                        cConv[key].append(dict())
                    cConv = cConv[key][index]
                else:
                    val = parts[-1]
                    if key.startswith(('al', 'an', 'ac')):
                        val = int(float(val))
                    elif key.startswith(('ad', 'afl')):
                        val = float(val)
                    elif key.startswith(('auc', 'aui', 'aul', 'aun')):
                        val = int(val, 16)
                    else:
                        val = val.strip('"')
                    for k in range(len(cConv[key]) - 1, index - 1):
                        cConv.get(key).append([])
                    cConv.get(key).insert(index, val)
            else:
                if i < len(splitparts) - 1:
                    if var not in cConv.keys():
                        cConv.update({var: dict()})
                    cConv = cConv[var]
                else:
                    key = splitparts[-1]
                    val = parts[-1]
                    if key.startswith(('l', 'n', 'c')):
                        val = int(float(val))
                    elif key.startswith(('d', 'fl')):
                        val = float(val)
                    elif key.startswith(('uc', 'ui', 'ul', 'un', 'e', 'size')):
                        val = int(val, 16)
                    else:
                        val = val.strip('"')
                    cConv.update({key: val})
    return ascconv


def parse_ascconv(string):
    """
    Parses xprotocol of siemens twix-file object and returns a protocol dict.
    """

    # output dict:
    ascconv = dict()
    for line in string.split('\n'):
        if line.startswith('#') or not line:
            continue
        cConv = ascconv
        parts = line.split()
        splitparts = parts[0].split('.')
        for i, var in enumerate(splitparts):
            if var.startswith('a'):
                tmp = var.strip(']').split('[')
                key = tmp[0]
                if len(tmp) > 1:
                    index = int(tmp[1])
                else:
                    index = 0
                if key not in cConv.keys():
                    cConv.update({key: list()})
                if i < len(splitparts) - 1:
                    for k in range(len(cConv[key]) - 1, index):
                        cConv[key].append(dict())
                    cConv = cConv[key][index]
                else:
                    val = parts[-1]
                    if key.startswith(('al', 'an', 'ac')):
                        val = int(float(val))
                    elif key.startswith(('ad', 'afl')):
                        val = float(val)
                    elif key.startswith(('auc', 'aui', 'aul', 'aun')):
                        val = int(val, 16)
                    else:
                        val = val.strip('"')
                    for k in range(len(cConv[key]) - 1, index - 1):
                        cConv.get(key).append([])
                    cConv.get(key).insert(index, val)
            else:
                if i < len(splitparts) - 1:
                    if var not in cConv.keys():
                        cConv.update({var: dict()})
                    cConv = cConv[var]
                else:
                    key = splitparts[-1]
                    val = parts[-1]
                    if key.startswith(('l', 'n', 'c')):
                        val = int(float(val))
                    elif key.startswith(('d', 'fl')):
                        val = float(val)
                    elif key.startswith(('uc', 'ui', 'ul', 'un', 'e', 'size')):
                        val = int(val, 16)
                    else:
                        val = val.strip('"')
                    cConv.update({key: val})
    return ascconv


def parse_xprot(buf):
    for line in buf.split('\n'):
        if not line.rsplit():
            continue
            
        
    return buf


def parse_twix_hdr(file):
    import re
    import struct
    import numpy as np
    hdr_len, n_buffer = np.fromfile(file, dtype=np.uint32, count=2)
    xprotocol = dict()
    pattern = b'(\w{4,})\x00(.{4})'
    pos = file.tell()
    for b in range(n_buffer):
        tmp = file.read(48)
        matches = re.search(pattern, tmp, re.DOTALL)
        name = matches.group(1).decode('latin1')
        buf_len = struct.unpack('<I', matches.group(2))[0]
        pos += len(matches.group())
        file.seek(pos)
        buf = file.read(buf_len).decode('latin1')
        if name == "MeasYaps":
            xprotocol[name] = parse_ascconv(buf)
        else:
            xprotocol[name] = parse_xprot(buf)
        pos += buf_len

    return xprotocol

    # header = fin.read(hdr_len)
    # header = header[:-24].decode('latin-1')
    # res = header.split("\n")
    # res = header
    # twixHeader = dict()
    # for linNum, line in enumerate(res):
    #     if line.lstrip().startswith("<ParamString"):
    #         name, value = readParamString(res, linNum, line)
    #         twixHeader[name] = value
    #     if line.lstrip().startswith("<ParamLong"):
    #         name, value = readParamLong(res, linNum, line)
    #         twixHeader[name] = value
    #     if line.lstrip().startswith("<ParamDouble"):
    #         name, value = readParamDouble(res, linNum, line)
    #         twixHeader[name] = value
    #     if line.lstrip().startswith("<ParamArray"):
    #         name, value = readParamArray(res, linNum, line)
    #         twixHeader[name] = value
    # 
    #     return twixHeader


def readParamDouble(res, linNum, line):
    name = (re.search("<.*>", line).group())
    name = name.split("\"")
    name = name[-2]

    value = (re.search("{.*}", line))
    if value is not None:
        try:
            value = re.split('{ |} |\s', value.group())
            value = float(value[-3])
        except:
            value = ""
    else:
        subLine = linNum+1
        value = ""
        tmp = res[subLine].lstrip()
        while "}" not in tmp:
            subLine = subLine + 1
            tmp = res[subLine].lstrip()
            if "<" not in tmp:
                value = value + tmp
        value = value[:-1]
        value = value.split()
        if len(value) == 0:
            value = 0
        elif len(value) == 1:
            value = float(value[0])
        else:
            value = np.array(value, dtype=np.float)

    return[name, value]


def readParamLong(res, linNum, line):
    name = (re.search("<.*>", line).group())
    name = name.split("\"")
    name = name[-2]

    value = (re.search("{.*}", line))
    if value is not None:
        try:
            value = re.split('{ |}', value.group())
            value = int(value[-2])
        except:
            value = ""
    else:
        subLine = linNum+1
        value = ""
        tmp = res[subLine].lstrip()
        while "}" not in tmp:
            subLine = subLine + 1
            tmp = res[subLine].lstrip()
            if "<" not in tmp:
                value = value + tmp
        value = value[:-1]
        value = value.split()
        if len(value) == 0:
            value = 0
        elif len(value) == 1:
            value = int(value[0])
        else:
            value = np.array(value, dtype=np.int)

    return[name, value]


def readParamString(res, linNum, line):
    name = (re.search("<.*>", line).group())
    name = name.split("\"")
    name = name[-2]

    value = (re.search("{.*}", line))
    if value is not None:
        value = value.group().split("\"")
        try:
            value = value[-2]
        except:
            value = ""
    else:
        subLine = linNum+1
        value = ''
        tmp = res[subLine].lstrip()
        while "}" not in tmp:
            subLine = subLine + 1
            tmp = res[subLine].lstrip()
            if "<" not in tmp:
                value = value + tmp
        value = value[:-1]
        value = value.split()
        if len(value) == 0:
            value = ''
        elif len(value) == 1:
            value = value[0]
        else:
            value = np.array(value)

    return [name, value]


def readParamArray(res, linNum, line):
    name = (re.search("<.*>", line).group())
    name = name.split("\"")
    name = name[-2]
    start = res[linNum+1]
    end = start.replace('{', '}')

    tmp = res[linNum+1]
    subLine = linNum + 1
    value = tmp.lstrip()
    while tmp != end:
        subLine = subLine + 1
        tmp = res[subLine]
        if tmp.lstrip() != '{ }':
            if "<" not in tmp:
                value = value + tmp.lstrip()

    value = re.sub("{|}", '', value)
    if value == '':
        value = 0

    return [name, value]
