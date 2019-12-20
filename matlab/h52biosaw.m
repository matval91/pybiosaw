function h52biosaw(coil)
    group = coil;
    x=h5read('efccdata.h5',['/' group '/coil_x']);
    y=h5read('efccdata.h5',['/' group '/coil_y']');
    z=h5read('efccdata.h5',['/' group '/coil_z']');
    str.X=x; str.Y=y; str.Z=z; 
    str.periodic=0;str.comments = {['JT60SA coil data. ' group]};
    addpath('/home/vallar/ASCOT/trunk/matlab/');
    write_ascot4_biosaw(str,['./' group 'coil_bio.txt']);
end