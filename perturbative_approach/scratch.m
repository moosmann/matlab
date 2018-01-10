>> whos
Name                         Size               Bytes  Class     Attributes

angles                       1x800               6400  double
gpu_index                    1x2                   16  double
link_data                    1x1                    8  double
pixel_size                   1x1                    8  double
rotation_axis_offset         1x1                    8  double
sino                      1280x800            4096000  single
tilt                         1x1                    8  double
tilt_lamino                  1x1                    8  double
vol_shape                    1x3                   24  double
vol_size                     1x6                   48  double


K>> sino(1)
single
-0.0383

K>> sino(end)
single
0.0266

K>> angles(1) = 0

K>> angles(end) = 6.2753

rotation_axis_offset = 580

vol_shape = 1280        1280           1

vol_size = -640.0000  640.0000 -640.0000  640.0000   -0.5000    0.5000

pixel_size = 1

link_data = 1

tilt = 0

gpu_index = 1     2

tilt_lamino = 0