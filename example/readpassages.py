import struct
import numpy as np

filename = 'passages1.in'

f = open(filename, 'rb')

headerlen = 12
a = f.read(headerlen)
header = struct.unpack('hhihh', a)

npasses = header[2]

linelen = 44
a = f.read(linelen*npasses)

lines = struct.unpack('2h9f2h'*npasses, a)
lines = np.asarray(lines)

inds = np.arange(0,npasses,dtype='int')
tenc = lines[inds*13 + 2]
time = lines[inds*13 + 3]
mass = lines[inds*13 + 4]
x = lines[inds*13 + 5]
y = lines[inds*13 + 6]
z = lines[inds*13 + 7]
vx = lines[inds*13 + 8]
vy = lines[inds*13 + 9]
vz = lines[inds*13 + 10]
rstart = np.sqrt(x**2.+y**2.+z**2.)

vmag = np.sqrt(vx**2.+vy**2. + vz**2.)
tbar = ((0.-x)*vx + (0.-y)*vy + (0.-z)*vz) / (vmag*vmag)
xclosebar = x + vx*tbar
yclosebar = y + vy*tbar
zclosebar = z + vz*tbar
rbar = np.sqrt(xclosebar**2. + yclosebar**2. + zclosebar**2.)

dsununitx = xclosebar/rbar
dsununity = yclosebar/rbar
dsununitz = zclosebar/rbar

impfactor = 2.*mass / vmag
impsunx = impfactor*dsununitx/rbar
impsuny = impfactor*dsununity/rbar
impsunz = impfactor*dsununitz/rbar
impsun = np.sqrt(impsunx**2.+impsuny**2.+impsunz**2.)
