# Fast racing line generation

Python3 implementation of the method described in https://arxiv.org/abs/1902.00606

+ It is written in <b>Python 3</b>

+  Python modules <b>required</b> are:
    + casadi
    + boostrtrees
    + scipy
    + json
    + numpy
    + matplotlib
    + os
    + sys
    + math
    + time
<br />

+ **boostrtrees** is downloadable from <a href="https://github.com/zouzias/pyboostrtree">here</a>. You can install it following these instructions:
    + Download the zip
    + Extract it
    + Install boost libraries. In Ubuntu 18.04 type `sudo apt-get install libboost-all-dev`. In ubuntu 16.04 you need to download and install it manually from <a href="https://sourceforge.net/projects/boost/files/latest/download">here</a> (it is needed a version >= 1.59)
    + Type in terminal `export BOOST_ROOT="/usr/include/boost"`. NOTE you need (prima devi aver installato boost)
    + Type `sudo pip3 install cython disttools` 
    + Inside the directory extracted from the zip type `python3 setup.py sdist`
    + Then type `pip3 install —upgrade path_libreria/dist/boostrtrees-0.0.1a1.tar.gz`
    <br />

+ The script **handles in input**:
    + `"path/<track_name>.json"` a string with the full path to the track's _.json_;
    + `"<track_name>"` a string containing the name of track; 
    + `<vMax>` max velocity (m/s);
    + `<gAcc>` max g, in module, in acceleration (g);
    + `<gDec` max g, in module, in deceleration (g);
    + `<gLat>` max lateral g (g);
    + `<safety_margin>` safety margin from the bounds (m);
    + `"<output_destination>"` a string containing the path where the output must be written.
    <br />

```
example of usage: 
python3 main.py ../SA_Modena_v2.1.json SA_Modena_v2.1 216 1.1 1.5 1.3 0.8 .

```
<br />

+ The script gives in output some files in the same directory:
    + `ft_gda_<track_name>.json`, containing racing line and speed profile
    + `racing_line<n_iteration>.png`, graphics containing the final racing line
