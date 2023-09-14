import mdtraj
import numpy as np

def get_rotcorrf(traj_prefix,trajectory_name, topology_name, vector_origin, vector_end, subtraj_length, dt):
    print("<<< Reading %s"%(trajectory_name))
    
    chunk_size = subtraj_length/dt    
   
    ## iterate trajectory in fragments
    for nchunks, chunk in enumerate(mdtraj.iterload(trajectory_name, top=topology_name,chunk=chunk_size)):
        chunk.unitcell_angles = None
        chunk.unitcell_lengths = None
        chunk.unitcell_vectors = None

        o_vec_idxs=chunk.topology.select(vector_origin)
        e_vec_idxs=chunk.topology.select(vector_end)
    
    
        unit_vec=np.zeros((len(o_vec_idxs),chunk.n_frames,3))
    
    #xyz is of shape (n_frames, n_atoms, 3)
        for i in range(len(chunk.xyz)):
            vec = chunk.xyz[i,e_vec_idxs,:] - chunk.xyz[i,o_vec_idxs,:]
            norm = np.atleast_2d(np.linalg.norm(vec, axis=1)).T
            unit_vec[:,i,:] = vec/norm
    
     
    
        Corr = np.zeros((len(o_vec_idxs),chunk.n_frames))
        print("<<< Calculating correlation function...")
        for i in range(0,len(o_vec_idxs)):
            for j in range(0,chunk.n_frames):
                # cosine of angle at different lag times is the scalar(dot) product of unit vectors
                # can be slow for large trajectories-> use numpy correlate instead
                scalar_product = np.sum(unit_vec[i,0:chunk.n_frames-j,:]*unit_vec[i,j:chunk.n_frames,:],axis=1) 
                Corr[i,j] = np.sum(0.5*(3*scalar_product**2-1),axis=0)/chunk.n_frames #Second Legendre polynomial
            
        AvgCorr=np.mean(Corr,axis=0) #average over all lysines within a chunk 
        np.savetxt(f'LYS_CB-NZ_rotcorrf_{traj_prefix}_'+str(nchunks)+'.txt',Corr.T, fmt='%.6f') #save individual Lys RotCorFun for each chunk
        np.savetxt(f'LYS_N-H_AVGrotcorrf_{traj_prefix}_'+str(nchunks)+'.txt',AvgCorr, fmt='%.6f')  #save average RotCorFun over all lysines within a chunk       
     



def main():
    for traj_prefix in ["prod1","prod2","prod3","prod4","prod5"]:
        trajectory_name = "traj_%s_not_aligned.dcd"%(traj_prefix) # for aligned traj, overall rotation of protein will be removed
        topology_name = "topology.pdb"
        
        vector_origin = "resname LYS and name N"
        vector_end = "resname LYS and name H"
        
        subtraj_length = 500 # in ns
        dt = 0.02 # in ns
            
        get_rotcorrf(traj_prefix,trajectory_name, topology_name, vector_origin, vector_end, subtraj_length, dt)
    
if __name__ == "__main__":
    main()

    