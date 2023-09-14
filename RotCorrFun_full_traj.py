import mdtraj
import numpy as np


def get_rotcorrf(traj_prefix,trajectory_name, topology_name, vector_origin, vector_end):
    print("<<< Reading %s"%(trajectory_name))
    
   
    ## load trajectory 
    traj = mdtraj.load(trajectory_name, top=topology_name)
    traj.unitcell_angles = None
    traj.unitcell_lengths = None
    traj.unitcell_vectors = None

    o_vec_idxs = traj.topology.select(vector_origin)
    e_vec_idxs = traj.topology.select(vector_end)
    
    
    unit_vec = np.zeros((len(o_vec_idxs),traj.n_frames,3))
    
    #xyz is of shape (n_frames, n_atoms, 3)
    for i in range(len(traj.xyz)):
        vec = traj.xyz[i,e_vec_idxs,:] - traj.xyz[i,o_vec_idxs,:]
        norm = np.atleast_2d(np.linalg.norm(vec,axis=1)).T
        unit_vec[:,i,:] = vec/norm
        
 
    
    Corr=np.zeros((len(o_vec_idxs),traj.n_frames)) 
    print("<<< Calculating correlation function...")
    for i in range(0,len(o_vec_idxs)):
        for j in range(0,traj.n_frames):
            # cosine of angle at different lag times is the scalar(dot) product of unit vectors
            # can be slow for large trajectories-> use numpy correlate instead
            scalar_product = np.sum(unit_vec[i,0:traj.n_frames-j,:]*unit_vec[i,j:traj.n_frames,:],axis=1)
            Corr[i,j] = np.sum(0.5*(3*scalar_product**2-1),axis=0)/traj.n_frames #Second Legendre polynomial
            
    AvgCorr=np.mean(Corr,axis=0) #average over all lysines
    np.savetxt(f'LYS_N-H_rotcorrf_{traj_prefix}.txt',Corr.T, fmt='%.6f') #save individual RotCorrFun for each of the lysines
    np.savetxt(f'LYS_N-H_AVGrotcorrf_{traj_prefix}.txt', AvgCorr, fmt='%.6f') #save average RotCorrFun over all lysines
   
     


def main():
    
   
    for traj_prefix in ["prod1","prod2","prod3","prod4","prod5"]:
        trajectory_name = "traj_%s_not_aligned.dcd"%(traj_prefix) # for aligned traj, overall rotation of protein will be removed
        topology_name = "topology.pdb"
        
        vector_origin = "resname LYS and name N"
        vector_end = "resname LYS and name H"
            
        get_rotcorrf(traj_prefix,trajectory_name, topology_name, vector_origin, vector_end)
    
if __name__ == "__main__":
    main()

    