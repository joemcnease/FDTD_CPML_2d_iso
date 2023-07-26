make
rm -d pressure/* vz/* vx/*
./main
python make_movie.py
