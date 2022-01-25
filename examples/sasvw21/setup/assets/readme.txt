How to create a ***vvfile***, which contains the lists of vessel (r,z)
Author: Alyssa Hayes and Ray Mattes

1) Ensure that this assets folder has a very important gfortran file, called gfile_vessel.f


2) Add the experiment G-file to this assets folder

mv /path/to/gfile/g123456.12345 /path/to/example/sasvw/run00/setup/assets


3) Create a symlink to the G-file

ln -s g123456.12345 gfile


4) Create an a.out file

gfortran gfile_vessel.f


5) Run a.out to get the vvfile

./a.out gfile
