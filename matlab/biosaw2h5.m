function biosaw2h5(folder, coil)
    coil=[coil 'coil'];
    addpath('/home/vallar/ASCOT/trunk/matlab/')
    f=read_ascot4([folder '/TFcoil.bpert']);
    
    outfile = [folder '/' coil 'out.h5'];
    h5create(outfile,['/Rgrid'],size(f.r(:,1)));
    h5write(outfile,['/Rgrid'],f.r);
    h5create(outfile,['/phigrid'],size(f.phi(:,1)));
    h5write(outfile,['/phigrid'],f.phi);    
    h5create(outfile,['/zgrid'],size(f.z(:,1)));
    h5write(outfile,['/zgrid'],f.z);
    
    h5create(outfile,['/' coil '/BR'],size(f.BR));
    h5write(outfile,['/' coil '/BR'],f.BR);    
    h5create(outfile,['/' coil '/Bphi'],size(f.Bphi));
    h5write(outfile,['/' coil '/Bphi'],f.Bphi);       
    h5create(outfile,['/' coil '/Bz'],size(f.Bz));
    h5write(outfile,['/' coil '/Bz'],f.Bz);       
    
    d=importdata([folder '/TFcoil_bio.txt'],' ',11);
    h5create(outfile,['/' coil '/coil_x'],size(d.data(:,1)));
    h5write(outfile,['/' coil '/coil_x'],d.data(:,1));
    h5create(outfile,['/' coil '/coil_y'],size(d.data(:,2)));
    h5write(outfile,['/' coil '/coil_y'],d.data(:,2));
    h5create(outfile,['/' coil '/coil_z'],size(d.data(:,3)));
    h5write(outfile,['/' coil '/coil_z'],d.data(:,3));
end