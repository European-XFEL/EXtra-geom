; JUNGFRAU geometry file written by EXtra-geom 1.2.1
; You may need to edit this file to add:
; - data and mask locations in the file
; - mask_good & mask_bad values to interpret the mask
; - adu_per_eV & photon_energy
; - clen (detector distance)
;
; See: http://www.desy.de/~twhite/crystfel/manual-crystfel_geometry.html

data = /data/data ;
dim0 = %
res = 13333.333333333334 ; pixels per metre

; Beam energy in eV
photon_energy = 4972

; Camera length, aka detector distance
clen = 0.101

; Analogue Digital Units per eV
adu_per_eV = 1
rigid_group_p0 = p0a0,p0a1,p0a2,p0a3,p0a4,p0a5,p0a6,p0a7
rigid_group_p1 = p1a0,p1a1,p1a2,p1a3,p1a4,p1a5,p1a6,p1a7
rigid_group_p2 = p2a0,p2a1,p2a2,p2a3,p2a4,p2a5,p2a6,p2a7
rigid_group_p3 = p3a0,p3a1,p3a2,p3a3,p3a4,p3a5,p3a6,p3a7
rigid_group_p4 = p4a0,p4a1,p4a2,p4a3,p4a4,p4a5,p4a6,p4a7
rigid_group_p5 = p5a0,p5a1,p5a2,p5a3,p5a4,p5a5,p5a6,p5a7
rigid_group_p6 = p6a0,p6a1,p6a2,p6a3,p6a4,p6a5,p6a6,p6a7
rigid_group_p7 = p7a0,p7a1,p7a2,p7a3,p7a4,p7a5,p7a6,p7a7
rigid_group_collection_modules = p0,p1,p2,p3,p4,p5,p6,p7

p0a0/dim1 = ss
p0a0/dim2 = fs
p0a0/min_fs = 0
p0a0/min_ss = 0
p0a0/max_fs = 255
p0a0/max_ss = 255
p0a0/fs = -1.0x +0.0y
p0a0/ss = +0.0x -1.0y
p0a0/corner_x = 1035.0
p0a0/corner_y = 1078.0
p0a0/coffset = 0.0

p0a1/dim1 = ss
p0a1/dim2 = fs
p0a1/min_fs = 256
p0a1/min_ss = 0
p0a1/max_fs = 511
p0a1/max_ss = 255
p0a1/fs = -1.0x +0.0y
p0a1/ss = +0.0x -1.0y
p0a1/corner_x = 777.0
p0a1/corner_y = 1078.0
p0a1/coffset = 0.0

p0a2/dim1 = ss
p0a2/dim2 = fs
p0a2/min_fs = 512
p0a2/min_ss = 0
p0a2/max_fs = 767
p0a2/max_ss = 255
p0a2/fs = -1.0x +0.0y
p0a2/ss = +0.0x -1.0y
p0a2/corner_x = 519.0
p0a2/corner_y = 1078.0
p0a2/coffset = 0.0

p0a3/dim1 = ss
p0a3/dim2 = fs
p0a3/min_fs = 768
p0a3/min_ss = 0
p0a3/max_fs = 1023
p0a3/max_ss = 255
p0a3/fs = -1.0x +0.0y
p0a3/ss = +0.0x -1.0y
p0a3/corner_x = 261.0
p0a3/corner_y = 1078.0
p0a3/coffset = 0.0

p0a4/dim1 = ss
p0a4/dim2 = fs
p0a4/min_fs = 0
p0a4/min_ss = 256
p0a4/max_fs = 255
p0a4/max_ss = 511
p0a4/fs = -1.0x +0.0y
p0a4/ss = +0.0x -1.0y
p0a4/corner_x = 1035.0
p0a4/corner_y = 820.0
p0a4/coffset = 0.0

p0a5/dim1 = ss
p0a5/dim2 = fs
p0a5/min_fs = 256
p0a5/min_ss = 256
p0a5/max_fs = 511
p0a5/max_ss = 511
p0a5/fs = -1.0x +0.0y
p0a5/ss = +0.0x -1.0y
p0a5/corner_x = 777.0
p0a5/corner_y = 820.0
p0a5/coffset = 0.0

p0a6/dim1 = ss
p0a6/dim2 = fs
p0a6/min_fs = 512
p0a6/min_ss = 256
p0a6/max_fs = 767
p0a6/max_ss = 511
p0a6/fs = -1.0x +0.0y
p0a6/ss = +0.0x -1.0y
p0a6/corner_x = 519.0
p0a6/corner_y = 820.0
p0a6/coffset = 0.0

p0a7/dim1 = ss
p0a7/dim2 = fs
p0a7/min_fs = 768
p0a7/min_ss = 256
p0a7/max_fs = 1023
p0a7/max_ss = 511
p0a7/fs = -1.0x +0.0y
p0a7/ss = +0.0x -1.0y
p0a7/corner_x = 261.0
p0a7/corner_y = 820.0
p0a7/coffset = 0.0

p1a0/dim1 = ss
p1a0/dim2 = fs
p1a0/min_fs = 0
p1a0/min_ss = 512
p1a0/max_fs = 255
p1a0/max_ss = 767
p1a0/fs = -1.0x +0.0y
p1a0/ss = +0.0x -1.0y
p1a0/corner_x = 1035.0
p1a0/corner_y = 531.0
p1a0/coffset = 0.0

p1a1/dim1 = ss
p1a1/dim2 = fs
p1a1/min_fs = 256
p1a1/min_ss = 512
p1a1/max_fs = 511
p1a1/max_ss = 767
p1a1/fs = -1.0x +0.0y
p1a1/ss = +0.0x -1.0y
p1a1/corner_x = 777.0
p1a1/corner_y = 531.0
p1a1/coffset = 0.0

p1a2/dim1 = ss
p1a2/dim2 = fs
p1a2/min_fs = 512
p1a2/min_ss = 512
p1a2/max_fs = 767
p1a2/max_ss = 767
p1a2/fs = -1.0x +0.0y
p1a2/ss = +0.0x -1.0y
p1a2/corner_x = 519.0
p1a2/corner_y = 531.0
p1a2/coffset = 0.0

p1a3/dim1 = ss
p1a3/dim2 = fs
p1a3/min_fs = 768
p1a3/min_ss = 512
p1a3/max_fs = 1023
p1a3/max_ss = 767
p1a3/fs = -1.0x +0.0y
p1a3/ss = +0.0x -1.0y
p1a3/corner_x = 261.0
p1a3/corner_y = 531.0
p1a3/coffset = 0.0

p1a4/dim1 = ss
p1a4/dim2 = fs
p1a4/min_fs = 0
p1a4/min_ss = 768
p1a4/max_fs = 255
p1a4/max_ss = 1023
p1a4/fs = -1.0x +0.0y
p1a4/ss = +0.0x -1.0y
p1a4/corner_x = 1035.0
p1a4/corner_y = 273.0
p1a4/coffset = 0.0

p1a5/dim1 = ss
p1a5/dim2 = fs
p1a5/min_fs = 256
p1a5/min_ss = 768
p1a5/max_fs = 511
p1a5/max_ss = 1023
p1a5/fs = -1.0x +0.0y
p1a5/ss = +0.0x -1.0y
p1a5/corner_x = 777.0
p1a5/corner_y = 273.0
p1a5/coffset = 0.0

p1a6/dim1 = ss
p1a6/dim2 = fs
p1a6/min_fs = 512
p1a6/min_ss = 768
p1a6/max_fs = 767
p1a6/max_ss = 1023
p1a6/fs = -1.0x +0.0y
p1a6/ss = +0.0x -1.0y
p1a6/corner_x = 519.0
p1a6/corner_y = 273.0
p1a6/coffset = 0.0

p1a7/dim1 = ss
p1a7/dim2 = fs
p1a7/min_fs = 768
p1a7/min_ss = 768
p1a7/max_fs = 1023
p1a7/max_ss = 1023
p1a7/fs = -1.0x +0.0y
p1a7/ss = +0.0x -1.0y
p1a7/corner_x = 261.0
p1a7/corner_y = 273.0
p1a7/coffset = 0.0

p2a0/dim1 = ss
p2a0/dim2 = fs
p2a0/min_fs = 0
p2a0/min_ss = 1024
p2a0/max_fs = 255
p2a0/max_ss = 1279
p2a0/fs = -1.0x +0.0y
p2a0/ss = +0.0x -1.0y
p2a0/corner_x = 1035.0
p2a0/corner_y = -16.0
p2a0/coffset = 0.0

p2a1/dim1 = ss
p2a1/dim2 = fs
p2a1/min_fs = 256
p2a1/min_ss = 1024
p2a1/max_fs = 511
p2a1/max_ss = 1279
p2a1/fs = -1.0x +0.0y
p2a1/ss = +0.0x -1.0y
p2a1/corner_x = 777.0
p2a1/corner_y = -16.0
p2a1/coffset = 0.0

p2a2/dim1 = ss
p2a2/dim2 = fs
p2a2/min_fs = 512
p2a2/min_ss = 1024
p2a2/max_fs = 767
p2a2/max_ss = 1279
p2a2/fs = -1.0x +0.0y
p2a2/ss = +0.0x -1.0y
p2a2/corner_x = 519.0
p2a2/corner_y = -16.0
p2a2/coffset = 0.0

p2a3/dim1 = ss
p2a3/dim2 = fs
p2a3/min_fs = 768
p2a3/min_ss = 1024
p2a3/max_fs = 1023
p2a3/max_ss = 1279
p2a3/fs = -1.0x +0.0y
p2a3/ss = +0.0x -1.0y
p2a3/corner_x = 261.0
p2a3/corner_y = -16.0
p2a3/coffset = 0.0

p2a4/dim1 = ss
p2a4/dim2 = fs
p2a4/min_fs = 0
p2a4/min_ss = 1280
p2a4/max_fs = 255
p2a4/max_ss = 1535
p2a4/fs = -1.0x +0.0y
p2a4/ss = +0.0x -1.0y
p2a4/corner_x = 1035.0
p2a4/corner_y = -274.0
p2a4/coffset = 0.0

p2a5/dim1 = ss
p2a5/dim2 = fs
p2a5/min_fs = 256
p2a5/min_ss = 1280
p2a5/max_fs = 511
p2a5/max_ss = 1535
p2a5/fs = -1.0x +0.0y
p2a5/ss = +0.0x -1.0y
p2a5/corner_x = 777.0
p2a5/corner_y = -274.0
p2a5/coffset = 0.0

p2a6/dim1 = ss
p2a6/dim2 = fs
p2a6/min_fs = 512
p2a6/min_ss = 1280
p2a6/max_fs = 767
p2a6/max_ss = 1535
p2a6/fs = -1.0x +0.0y
p2a6/ss = +0.0x -1.0y
p2a6/corner_x = 519.0
p2a6/corner_y = -274.0
p2a6/coffset = 0.0

p2a7/dim1 = ss
p2a7/dim2 = fs
p2a7/min_fs = 768
p2a7/min_ss = 1280
p2a7/max_fs = 1023
p2a7/max_ss = 1535
p2a7/fs = -1.0x +0.0y
p2a7/ss = +0.0x -1.0y
p2a7/corner_x = 261.0
p2a7/corner_y = -274.0
p2a7/coffset = 0.0

p3a0/dim1 = ss
p3a0/dim2 = fs
p3a0/min_fs = 0
p3a0/min_ss = 1536
p3a0/max_fs = 255
p3a0/max_ss = 1791
p3a0/fs = -1.0x +0.0y
p3a0/ss = +0.0x -1.0y
p3a0/corner_x = 1035.0
p3a0/corner_y = -563.0
p3a0/coffset = 0.0

p3a1/dim1 = ss
p3a1/dim2 = fs
p3a1/min_fs = 256
p3a1/min_ss = 1536
p3a1/max_fs = 511
p3a1/max_ss = 1791
p3a1/fs = -1.0x +0.0y
p3a1/ss = +0.0x -1.0y
p3a1/corner_x = 777.0
p3a1/corner_y = -563.0
p3a1/coffset = 0.0

p3a2/dim1 = ss
p3a2/dim2 = fs
p3a2/min_fs = 512
p3a2/min_ss = 1536
p3a2/max_fs = 767
p3a2/max_ss = 1791
p3a2/fs = -1.0x +0.0y
p3a2/ss = +0.0x -1.0y
p3a2/corner_x = 519.0
p3a2/corner_y = -563.0
p3a2/coffset = 0.0

p3a3/dim1 = ss
p3a3/dim2 = fs
p3a3/min_fs = 768
p3a3/min_ss = 1536
p3a3/max_fs = 1023
p3a3/max_ss = 1791
p3a3/fs = -1.0x +0.0y
p3a3/ss = +0.0x -1.0y
p3a3/corner_x = 261.0
p3a3/corner_y = -563.0
p3a3/coffset = 0.0

p3a4/dim1 = ss
p3a4/dim2 = fs
p3a4/min_fs = 0
p3a4/min_ss = 1792
p3a4/max_fs = 255
p3a4/max_ss = 2047
p3a4/fs = -1.0x +0.0y
p3a4/ss = +0.0x -1.0y
p3a4/corner_x = 1035.0
p3a4/corner_y = -821.0
p3a4/coffset = 0.0

p3a5/dim1 = ss
p3a5/dim2 = fs
p3a5/min_fs = 256
p3a5/min_ss = 1792
p3a5/max_fs = 511
p3a5/max_ss = 2047
p3a5/fs = -1.0x +0.0y
p3a5/ss = +0.0x -1.0y
p3a5/corner_x = 777.0
p3a5/corner_y = -821.0
p3a5/coffset = 0.0

p3a6/dim1 = ss
p3a6/dim2 = fs
p3a6/min_fs = 512
p3a6/min_ss = 1792
p3a6/max_fs = 767
p3a6/max_ss = 2047
p3a6/fs = -1.0x +0.0y
p3a6/ss = +0.0x -1.0y
p3a6/corner_x = 519.0
p3a6/corner_y = -821.0
p3a6/coffset = 0.0

p3a7/dim1 = ss
p3a7/dim2 = fs
p3a7/min_fs = 768
p3a7/min_ss = 1792
p3a7/max_fs = 1023
p3a7/max_ss = 2047
p3a7/fs = -1.0x +0.0y
p3a7/ss = +0.0x -1.0y
p3a7/corner_x = 261.0
p3a7/corner_y = -821.0
p3a7/coffset = 0.0

p4a0/dim1 = ss
p4a0/dim2 = fs
p4a0/min_fs = 0
p4a0/min_ss = 2048
p4a0/max_fs = 255
p4a0/max_ss = 2303
p4a0/fs = +1.0x +0.0y
p4a0/ss = +0.0x +1.0y
p4a0/corner_x = -1035.0
p4a0/corner_y = -1078.0
p4a0/coffset = 0.0

p4a1/dim1 = ss
p4a1/dim2 = fs
p4a1/min_fs = 256
p4a1/min_ss = 2048
p4a1/max_fs = 511
p4a1/max_ss = 2303
p4a1/fs = +1.0x +0.0y
p4a1/ss = +0.0x +1.0y
p4a1/corner_x = -777.0
p4a1/corner_y = -1078.0
p4a1/coffset = 0.0

p4a2/dim1 = ss
p4a2/dim2 = fs
p4a2/min_fs = 512
p4a2/min_ss = 2048
p4a2/max_fs = 767
p4a2/max_ss = 2303
p4a2/fs = +1.0x +0.0y
p4a2/ss = +0.0x +1.0y
p4a2/corner_x = -519.0
p4a2/corner_y = -1078.0
p4a2/coffset = 0.0

p4a3/dim1 = ss
p4a3/dim2 = fs
p4a3/min_fs = 768
p4a3/min_ss = 2048
p4a3/max_fs = 1023
p4a3/max_ss = 2303
p4a3/fs = +1.0x +0.0y
p4a3/ss = +0.0x +1.0y
p4a3/corner_x = -261.0
p4a3/corner_y = -1078.0
p4a3/coffset = 0.0

p4a4/dim1 = ss
p4a4/dim2 = fs
p4a4/min_fs = 0
p4a4/min_ss = 2304
p4a4/max_fs = 255
p4a4/max_ss = 2559
p4a4/fs = +1.0x +0.0y
p4a4/ss = +0.0x +1.0y
p4a4/corner_x = -1035.0
p4a4/corner_y = -820.0
p4a4/coffset = 0.0

p4a5/dim1 = ss
p4a5/dim2 = fs
p4a5/min_fs = 256
p4a5/min_ss = 2304
p4a5/max_fs = 511
p4a5/max_ss = 2559
p4a5/fs = +1.0x +0.0y
p4a5/ss = +0.0x +1.0y
p4a5/corner_x = -777.0
p4a5/corner_y = -820.0
p4a5/coffset = 0.0

p4a6/dim1 = ss
p4a6/dim2 = fs
p4a6/min_fs = 512
p4a6/min_ss = 2304
p4a6/max_fs = 767
p4a6/max_ss = 2559
p4a6/fs = +1.0x +0.0y
p4a6/ss = +0.0x +1.0y
p4a6/corner_x = -519.0
p4a6/corner_y = -820.0
p4a6/coffset = 0.0

p4a7/dim1 = ss
p4a7/dim2 = fs
p4a7/min_fs = 768
p4a7/min_ss = 2304
p4a7/max_fs = 1023
p4a7/max_ss = 2559
p4a7/fs = +1.0x +0.0y
p4a7/ss = +0.0x +1.0y
p4a7/corner_x = -261.0
p4a7/corner_y = -820.0
p4a7/coffset = 0.0

p5a0/dim1 = ss
p5a0/dim2 = fs
p5a0/min_fs = 0
p5a0/min_ss = 2560
p5a0/max_fs = 255
p5a0/max_ss = 2815
p5a0/fs = +1.0x +0.0y
p5a0/ss = +0.0x +1.0y
p5a0/corner_x = -1035.0
p5a0/corner_y = -531.0
p5a0/coffset = 0.0

p5a1/dim1 = ss
p5a1/dim2 = fs
p5a1/min_fs = 256
p5a1/min_ss = 2560
p5a1/max_fs = 511
p5a1/max_ss = 2815
p5a1/fs = +1.0x +0.0y
p5a1/ss = +0.0x +1.0y
p5a1/corner_x = -777.0
p5a1/corner_y = -531.0
p5a1/coffset = 0.0

p5a2/dim1 = ss
p5a2/dim2 = fs
p5a2/min_fs = 512
p5a2/min_ss = 2560
p5a2/max_fs = 767
p5a2/max_ss = 2815
p5a2/fs = +1.0x +0.0y
p5a2/ss = +0.0x +1.0y
p5a2/corner_x = -519.0
p5a2/corner_y = -531.0
p5a2/coffset = 0.0

p5a3/dim1 = ss
p5a3/dim2 = fs
p5a3/min_fs = 768
p5a3/min_ss = 2560
p5a3/max_fs = 1023
p5a3/max_ss = 2815
p5a3/fs = +1.0x +0.0y
p5a3/ss = +0.0x +1.0y
p5a3/corner_x = -261.0
p5a3/corner_y = -531.0
p5a3/coffset = 0.0

p5a4/dim1 = ss
p5a4/dim2 = fs
p5a4/min_fs = 0
p5a4/min_ss = 2816
p5a4/max_fs = 255
p5a4/max_ss = 3071
p5a4/fs = +1.0x +0.0y
p5a4/ss = +0.0x +1.0y
p5a4/corner_x = -1035.0
p5a4/corner_y = -273.0
p5a4/coffset = 0.0

p5a5/dim1 = ss
p5a5/dim2 = fs
p5a5/min_fs = 256
p5a5/min_ss = 2816
p5a5/max_fs = 511
p5a5/max_ss = 3071
p5a5/fs = +1.0x +0.0y
p5a5/ss = +0.0x +1.0y
p5a5/corner_x = -777.0
p5a5/corner_y = -273.0
p5a5/coffset = 0.0

p5a6/dim1 = ss
p5a6/dim2 = fs
p5a6/min_fs = 512
p5a6/min_ss = 2816
p5a6/max_fs = 767
p5a6/max_ss = 3071
p5a6/fs = +1.0x +0.0y
p5a6/ss = +0.0x +1.0y
p5a6/corner_x = -519.0
p5a6/corner_y = -273.0
p5a6/coffset = 0.0

p5a7/dim1 = ss
p5a7/dim2 = fs
p5a7/min_fs = 768
p5a7/min_ss = 2816
p5a7/max_fs = 1023
p5a7/max_ss = 3071
p5a7/fs = +1.0x +0.0y
p5a7/ss = +0.0x +1.0y
p5a7/corner_x = -261.0
p5a7/corner_y = -273.0
p5a7/coffset = 0.0

p6a0/dim1 = ss
p6a0/dim2 = fs
p6a0/min_fs = 0
p6a0/min_ss = 3072
p6a0/max_fs = 255
p6a0/max_ss = 3327
p6a0/fs = +1.0x +0.0y
p6a0/ss = +0.0x +1.0y
p6a0/corner_x = -1035.0
p6a0/corner_y = 16.0
p6a0/coffset = 0.0

p6a1/dim1 = ss
p6a1/dim2 = fs
p6a1/min_fs = 256
p6a1/min_ss = 3072
p6a1/max_fs = 511
p6a1/max_ss = 3327
p6a1/fs = +1.0x +0.0y
p6a1/ss = +0.0x +1.0y
p6a1/corner_x = -777.0
p6a1/corner_y = 16.0
p6a1/coffset = 0.0

p6a2/dim1 = ss
p6a2/dim2 = fs
p6a2/min_fs = 512
p6a2/min_ss = 3072
p6a2/max_fs = 767
p6a2/max_ss = 3327
p6a2/fs = +1.0x +0.0y
p6a2/ss = +0.0x +1.0y
p6a2/corner_x = -519.0
p6a2/corner_y = 16.0
p6a2/coffset = 0.0

p6a3/dim1 = ss
p6a3/dim2 = fs
p6a3/min_fs = 768
p6a3/min_ss = 3072
p6a3/max_fs = 1023
p6a3/max_ss = 3327
p6a3/fs = +1.0x +0.0y
p6a3/ss = +0.0x +1.0y
p6a3/corner_x = -261.0
p6a3/corner_y = 16.0
p6a3/coffset = 0.0

p6a4/dim1 = ss
p6a4/dim2 = fs
p6a4/min_fs = 0
p6a4/min_ss = 3328
p6a4/max_fs = 255
p6a4/max_ss = 3583
p6a4/fs = +1.0x +0.0y
p6a4/ss = +0.0x +1.0y
p6a4/corner_x = -1035.0
p6a4/corner_y = 274.0
p6a4/coffset = 0.0

p6a5/dim1 = ss
p6a5/dim2 = fs
p6a5/min_fs = 256
p6a5/min_ss = 3328
p6a5/max_fs = 511
p6a5/max_ss = 3583
p6a5/fs = +1.0x +0.0y
p6a5/ss = +0.0x +1.0y
p6a5/corner_x = -777.0
p6a5/corner_y = 274.0
p6a5/coffset = 0.0

p6a6/dim1 = ss
p6a6/dim2 = fs
p6a6/min_fs = 512
p6a6/min_ss = 3328
p6a6/max_fs = 767
p6a6/max_ss = 3583
p6a6/fs = +1.0x +0.0y
p6a6/ss = +0.0x +1.0y
p6a6/corner_x = -519.0
p6a6/corner_y = 274.0
p6a6/coffset = 0.0

p6a7/dim1 = ss
p6a7/dim2 = fs
p6a7/min_fs = 768
p6a7/min_ss = 3328
p6a7/max_fs = 1023
p6a7/max_ss = 3583
p6a7/fs = +1.0x +0.0y
p6a7/ss = +0.0x +1.0y
p6a7/corner_x = -261.0
p6a7/corner_y = 274.0
p6a7/coffset = 0.0

p7a0/dim1 = ss
p7a0/dim2 = fs
p7a0/min_fs = 0
p7a0/min_ss = 3584
p7a0/max_fs = 255
p7a0/max_ss = 3839
p7a0/fs = +1.0x +0.0y
p7a0/ss = +0.0x +1.0y
p7a0/corner_x = -1035.0
p7a0/corner_y = 563.0
p7a0/coffset = 0.0

p7a1/dim1 = ss
p7a1/dim2 = fs
p7a1/min_fs = 256
p7a1/min_ss = 3584
p7a1/max_fs = 511
p7a1/max_ss = 3839
p7a1/fs = +1.0x +0.0y
p7a1/ss = +0.0x +1.0y
p7a1/corner_x = -777.0
p7a1/corner_y = 563.0
p7a1/coffset = 0.0

p7a2/dim1 = ss
p7a2/dim2 = fs
p7a2/min_fs = 512
p7a2/min_ss = 3584
p7a2/max_fs = 767
p7a2/max_ss = 3839
p7a2/fs = +1.0x +0.0y
p7a2/ss = +0.0x +1.0y
p7a2/corner_x = -519.0
p7a2/corner_y = 563.0
p7a2/coffset = 0.0

p7a3/dim1 = ss
p7a3/dim2 = fs
p7a3/min_fs = 768
p7a3/min_ss = 3584
p7a3/max_fs = 1023
p7a3/max_ss = 3839
p7a3/fs = +1.0x +0.0y
p7a3/ss = +0.0x +1.0y
p7a3/corner_x = -261.0
p7a3/corner_y = 563.0
p7a3/coffset = 0.0

p7a4/dim1 = ss
p7a4/dim2 = fs
p7a4/min_fs = 0
p7a4/min_ss = 3840
p7a4/max_fs = 255
p7a4/max_ss = 4095
p7a4/fs = +1.0x +0.0y
p7a4/ss = +0.0x +1.0y
p7a4/corner_x = -1035.0
p7a4/corner_y = 821.0
p7a4/coffset = 0.0

p7a5/dim1 = ss
p7a5/dim2 = fs
p7a5/min_fs = 256
p7a5/min_ss = 3840
p7a5/max_fs = 511
p7a5/max_ss = 4095
p7a5/fs = +1.0x +0.0y
p7a5/ss = +0.0x +1.0y
p7a5/corner_x = -777.0
p7a5/corner_y = 821.0
p7a5/coffset = 0.0

p7a6/dim1 = ss
p7a6/dim2 = fs
p7a6/min_fs = 512
p7a6/min_ss = 3840
p7a6/max_fs = 767
p7a6/max_ss = 4095
p7a6/fs = +1.0x +0.0y
p7a6/ss = +0.0x +1.0y
p7a6/corner_x = -519.0
p7a6/corner_y = 821.0
p7a6/coffset = 0.0

p7a7/dim1 = ss
p7a7/dim2 = fs
p7a7/min_fs = 768
p7a7/min_ss = 3840
p7a7/max_fs = 1023
p7a7/max_ss = 4095
p7a7/fs = +1.0x +0.0y
p7a7/ss = +0.0x +1.0y
p7a7/corner_x = -261.0
p7a7/corner_y = 821.0
p7a7/coffset = 0.0
