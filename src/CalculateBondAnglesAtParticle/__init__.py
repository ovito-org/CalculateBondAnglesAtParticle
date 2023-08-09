### Calculate Bond Angles At Particle
# This modifier outputs the angles between all pairwise combinations of bonds at one particle.

from ovito.data import *
import numpy as np 
from traits.api import Union, Int, Enum, Bool
from ovito.pipeline import ModifierInterface
from itertools import combinations

class CalculateBondAnglesAtParticle(ModifierInterface):
   
    center_particle = Union(None, Int, label = "Compute for particle")
    mode = Enum("Index", "Identifier", label = "Choose particle by")
    bond_mode = Enum("Index", "Identifier", label = "List bonds by")

    def calculate_bond_angles(self, bond_vectors):
        b1 = bond_vectors[:,0]/np.linalg.norm(bond_vectors[:,0], axis = 1)[:, None]
        b2 = bond_vectors[:,1]/np.linalg.norm(bond_vectors[:,1], axis = 1)[:, None]
        return np.degrees(np.arccos(np.sum(b1*b2, axis = 1)))
    
    def calculate_bond_vector_combinations(self, data, particle):
        # Create bonds enumerator object.
        bonds_enum = BondsEnumerator(data.particles.bonds)
        # List of bond indices of all bonds at the current particle
        bonds_of_particle = [bond_index for bond_index in bonds_enum.bonds_of_particle(particle)]
        if len(bonds_of_particle) < 2:
            raise RuntimeError("Not enough bonds found to compute angles.")
        # All possible bond pairs in bonds_of_particle list
        bond_pairs = list(combinations(bonds_of_particle, 2))
        # Look up corresponding bond vectors
        topo = data.particles.bonds.topology[bond_pairs] 
        # Flip bond vector if current particle is not index 0 in topology pairs
        idx = np.where(topo[:,:,0] != particle)
        vectors = data.particles.bonds["Bond vectors"][bond_pairs]
        vectors[idx[0], idx[1]] *= -1
        # Get particle index triplets from topo [B,A][B,C]
        topo[idx[0], idx[1]] =  np.flip(topo[idx[0], idx[1]])
        triplets = np.column_stack((topo[:, 0, 1],topo[:, 1, :]))
        return bond_pairs, triplets, vectors
       
    def modify(self, data, **kwargs):

        if self.center_particle == None:
            return
        if self.mode == "Identifier" and "Particle Identifier" not in data.particles:
            raise RuntimeError("No Particle Identifiers in DataCollection. Deactivate Option <Output Particle Identifiers>.")
        if self.mode == "Identifier" and (np.isin(data.particles.identifiers, self.center_particle).any() == False):
            raise IndexError("Invalid Particle Identifier")
        if self.mode == "Index" and (self.center_particle >= data.particles.count or self.center_particle < 0):
            raise IndexError(f"Invalid Particle Index. Choose Index between 0 and {data.particles.count-1}.")
        if data.particles.bonds == None:
            raise RuntimeError("No Bonds in DataCollection. Please first generate bonds.")
        if self.bond_mode == "Identifier" and "Bond Identifier" not in data.particles.bonds:
            raise RuntimeError("No Bond Identifiers in DataCollection. Switch to Option <Output Bond Indices>.")

        # Look up bond indices of all bonds connected to current particle
        if self.mode == 'Identifier':
            particle = np.where(data.particles.identifiers == self.center_particle)[0][0]
        else:
            particle = self.center_particle

        # Calculate bond vectors for global topology array    
        positions = data.particles.positions
        topology = data.particles.bonds.topology
        bond_vectors = positions[topology[:,1]] - positions[topology[:,0]]
        if "Periodic Image" in data.particles.bonds:
            bond_vectors += np.dot(data.cell[:3,:3], data.particles.bonds.pbc_vectors.T).T
        data.particles_.bonds_.create_property("Bond vectors", data = bond_vectors)    
        # Get all possible combinations of bond pairs at one particle, the corresponding 
        # particle index triplets and bond vectors
        bond_pairs, triplets, v_b = self.calculate_bond_vector_combinations(data, particle)
        # Calculcate angles between all pairs of bond vectors
        angles = self.calculate_bond_angles(v_b)
        
        #Script output
        b = len(str(data.particles.bonds.count))
        p = len(str(data.particles.count))
        header = f"Triplet A-B-C{(2*p+10+len(str(particle))-13)*' '}Angle{7*' '}Bond Pair B1 - B2" 
        print(header)
        print(f"{'-'*len(header)}")
        for i in range(len(bond_pairs)):
            A,B,C = triplets[i]
            if self.mode == 'Identifier':
                A,B,C = data.particles.identifiers[[A, B, C]]
            output = f"{A:{p}d} - {B} - {C:<{p}d}    {angles[i]:<8.4f}" 
            B1, B2 = bond_pairs[i][0], bond_pairs[i][1]
            if self.bond_mode == 'Identifier':
                B1, B2 = int(data.particles.bonds['Bond Identifier'][B1]),int(data.particles.bonds['Bond Identifier'][B2])
            output+=f"    {B1:{b}d} - {B2:<{b}d}"
            print(output)
        
        # Store results as OVITO Data Table
        table = data.tables.create(
            identifier=f"bond-angles-{self.center_particle}",
            title=f"Bond Angles of Particle {self.center_particle}",
            plot_mode=DataTable.PlotMode.NoPlot)
        table.y = table.create_property('Angle', data=angles)
        if self.mode == 'Identifier':
            triplets = data.particles.identifiers[triplets]
        table.create_property('Particle Triplet', data=triplets, components=['A', 'B', 'C'])
    
        if self.bond_mode == 'Identifier':
            bond_pairs = data.particles.bonds['Bond Identifier'][bond_pairs]
        table.create_property('Bond Pair', data=bond_pairs, components=['Bond1', 'Bond2'])
