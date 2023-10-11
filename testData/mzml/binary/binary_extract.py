import base64
import struct
import re

file = "binary_mz"

bin = open(file,'rb')
encoded_data = bin.read()

## have to delete the training /n
#re.sub("/n$", "", encoded_data
encoded_data = encoded_data[:-1]

raw_data = zlib.decompress(base64.decodebytes(encoded_data))
out64_mz = struct.unpack("<%sd" % (len(raw_data) // 8), raw_data)
out32 = struct.unpack('<%sf' % (len(raw_data) // 4), raw_data)
