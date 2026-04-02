
# The code below can be run in R to perform the camera calibration, 
# After the calibration videos (left.mp4 and right.mp4) have been copied to the "calibration" directory, 
# and the observation videos (left.mp4 and right.mp4) have been copied to the "observations" directory.

calibrateCameras(img.dir = cal_dir, 
                 nx = 8,               # number of horizontal inner corners 
                 ny = 6,               # number of vertical inner corners  
                 sq.size = "3.04 cm", # square size
                 cal.file = 'calibration_gopro.txt',
                 corner.dir = 'corners', verify.dir = 'verify', error.dir = 'errors',
                 undistort = TRUE, num.aspects.read = 80, fit.min.break = 2, nlm.calls.max = 15,
                 objective.min = 0.8, max.sample.optim = 30, num.sample.est = 20,
                 num.aspects.sample = 4, num.sample.sets = 3, objective.min.break = 1.2)


# The code below can be run in R to digitze the images, 
# after the observation videos (left.mp4 and right.mp4) have been copied to the observation directory. 

extractFramesDir(img.dir = obs_dir, nth = 1)
file.copy("calibration_gopro.txt" , file.path(obs_dir,"calibration_gopro.txt"))

digitizeImages(image.file = file.path(obs_dir, 'Images'), 
               shapes.file = file.path(obs_dir, 'Shapes 2D'),
               landmarks.ref = file.path(obs_dir, 'landmarks.txt'), 
               cal.file = file.path(obs_dir,'calibration_gopro.txt'))


# The code below can be run in R to calculate the 3D coordinates of the landmarks, and the distance between 
# the snout and caudal fin landmarks. 

shapes <- reconstructStereoSets(shapes.2d = file.path(obs_dir, 'Shapes 2D'), 
                      shapes.3d = file.path(obs_dir, 'Shapes 3D'),
                      cal.file  = file.path(obs_dir, 'calibration_gopro.txt'))


distances(shapes)
