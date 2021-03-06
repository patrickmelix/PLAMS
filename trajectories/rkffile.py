#!/usr/bin/env python

import numpy

from scm.plams import KFFile, Molecule, Bond, AMSResults, Units, PeriodicTable, PlamsError
from .trajectoryfile import TrajectoryFile

__all__ = ['RKFTrajectoryFile']

bohr_to_angstrom = Units.conversion_ratio('bohr','angstrom')

class RKFTrajectoryFile (TrajectoryFile) :
        """
        Class that represents an RKF file
        """
        def __init__ (self, filename, mode='rb', fileobject=None, ntap=None) :
                """
                Creates an RKF file object
                """
                self.position = 0
                if filename is not None :
                        fileobject = KFFile(filename,autosave=False)
                        #fileobject = KFFile(filename,autosave=False,fastsave=True)
                        if fileobject is None :
                                raise PlamsError('KFFile %s not found.'%(rkfname))
                self.file_object = fileobject
                self.mode = mode

                self.ntap = 0
                if ntap is not None :
                        self.ntap = ntap
                self.firsttime = True
                self.coords = numpy.zeros((self.ntap,3))                # Only for reading purposes, to avoid creating the array each time

                # RKF specific attributes
                #self.saving_freq = 100
                self.nvecs = 3
                self.latticevecs = numpy.zeros((3,3))
                self.elements = ['H']*self.ntap
                self.read_lattice = True               # Reading time can be saved by skipping the lattice info
                self.read_bonds = True
                self.cell = numpy.zeros((3,3))
                self.conect = None
                self.timestep = None
                self.include_mddata = False
                self.mddata = None

                # Skip to the trajectory part of the file (only if in read mode, because coords are required in header)
                if self.mode == 'rb' :
                        self.read_header()
                elif self.mode == 'wb' :
                        sections = self.file_object.sections()
                        if len(sections) > 0 : 
                                raise PlamsError ('RKF file %s already exists'%(filename))
                else :
                        raise PlamsError ('Mode %s is invalid. Only "rb" and "wb" are allowed.'%(self.mode))

        def get_elements (self) :
                """
                Get the elements attribute
                """ 
                return self.elements

        def set_elements (self, elements) :
                """
                Sets the elements. Needed in write mode.

                @param elements : The element names of the atoms
                """
                self.elements = elements

        def set_name (self, name) :
                """
                Sets the name of the system. Needed in write mode
                """
                self.name = name

        def close (self) :
                """
                Makes sure that all commands that come before are executed.
                """
                # Write the step info
                if self.timestep is not None and self.mode == 'wb' :
                        self.file_object.write('MDResults','StartStep',0)
                        self.file_object.write('MDResults','StartTime[fs]',0.)
                        nsteps = self.get_length()
                        self.file_object.write('MDResults','EndStep',nsteps-1)
                        self.file_object.write('MDResults','EndTime[fs]',(nsteps-1)*self.timestep)

                # Write to file
                if self.mode == 'wb' :
                        self.file_object.save()
                del(self)

        def read_header (self) :
                """
                Set up info required for reading frames
                """
                self.elements = self.file_object.read('Molecule','AtomSymbols').split()
                self.elements = [el.split('.')[0] for el in self.elements]
                if 'MDHistory' in self.file_object.reader._sections :
                        self.timestep = self.file_object.read('MDHistory','Time(1)')[1]
                self.ntap = len(self.elements)
                self.coords = numpy.zeros((self.ntap,3))
                try :
                        self.latticevecs = numpy.array(self.file_object.read('Molecule','LatticeVectors'))
                        #self.nvecs = int(len(self.cell)/3)
                        # New code 27-05-2020
                        self.nvecs = int(len(self.latticevecs)/3)
                        self.latticevecs = self.latticevecs.reshape((self.nvecs,3))
                except KeyError :
                        pass

        def write_header (self, coords, cell) :
                """
                Write Molecule info to file (elements, periodicity)
                """
                # First write the general section
                self.file_object.write('General','file-ident','RKF')
                self.file_object.write('General','termination status','NORMAL TERMINATION')
                self.file_object.write('General','program','ams')

                # Then write the input molecule
                self._write_molecule_section(coords, cell)
                if self.include_mddata :
                        self._write_molecule_section(coords, cell, section='InputMolecule')
                        # Start setting up the MDHistory section as well
                        self.file_object.write('MDHistory','blockSize',100)

        def _write_molecule_section (self, coords, cell, section='Molecule') :
                """
                Write the molecule section
                """
                # Then write the input molecule
                charge = 0.
                element_numbers = [PeriodicTable.get_atomic_number(el) for el in self.elements]

                self.file_object.write(section,'nAtoms',self.ntap)
                self.file_object.write(section,'AtomicNumbers',element_numbers)
                self.file_object.write(section,'AtomSymbols',self.elements)
                crd = [Units.convert(float(c),'angstrom','bohr') for coord in coords for c in coord]
                self.file_object.write(section,'Coords',crd)
                self.file_object.write(section,'Charge',charge)
                if cell is not None :
                        self.file_object.write(section,'nLatticeVectors',self.nvecs)
                        vecs = [Units.convert(float(v),'angstrom','bohr') for vec in cell for v in vec]
                        self.file_object.write(section,'LatticeVectors',vecs)
                # Should it write bonds?

        def get_plamsmol (self) :
                """
                Creates a PLAMS molecule object from the xyz-trajectory file
                """
                section_dict = self.file_object.read_section('Molecule')
                plamsmol = Molecule._mol_from_rkf_section(section_dict)
                return plamsmol

        def read_frame (self, i, molecule=None) :
                """
                Reads the relevant info from frame i

                @param latticevecs: Numpy array for the cellvectors in kffile (self.nvecs*3)
                @param vecs       : Numpy array for output cellvectors (3,3)
                """
                # Read the coordinates
                try :
                        if not self.coords.shape == (self.ntap,3) :
                                raise PlamsError('coords attribute has been changed outside the class')
                        coords = self.coords.reshape(self.ntap*3)
                        coords[:] = self.file_object.read('History', 'Coords(%i)'%(i+1))
                        coords *= bohr_to_angstrom
                        # This has changed self.coords behind the scenes
                except KeyError :
                        return None, None
                
                # Read the cell data
                cell = None
                if self.read_lattice :
                        try :
                                latticevecs = self.latticevecs.reshape(self.nvecs*3)
                                latticevecs[:] = self.file_object.read('History','LatticeVectors(%i)'%(i+1)) #* bohr_to_angstrom
                                latticevecs *= bohr_to_angstrom
                                # This changed self.latticevecs behind the scenes
                                #self.cell[:self.nvecs] = latticevecs
                                self.cell[:self.nvecs] = self.latticevecs
                                cell = self.cell
                        except KeyError :
                                pass

                # Read the bond data
                conect = None
                bonds = None
                if self.read_bonds :
                        conect, bonds = self.read_bond_data(section='History', step=i)
                self.conect = conect

                # Read and store all MDData for this frame
                if self.include_mddata :
                        self.store_mddata_for_step(i)

                if isinstance(molecule,Molecule) :
                        self._set_plamsmol(self.coords,cell,molecule,bonds)                
               
                self.position = i
                return self.coords, cell

        def read_bond_data (self, section, step=None) :
                """
                Read the bond data from the rkf file
                """
                conect = None
                bonds = None
                try :
                        step_txt = ''
                        if step is not None :
                                step_txt = '(%i)'%(step+1)
                        indices = self.file_object.read(section,'Bonds.Index%s'%(step_txt))
                        connection_table = self.file_object.read(section,'Bonds.Atoms%s'%(step_txt))
                        if isinstance(connection_table,int) :
                                connection_table = [connection_table]
                        bond_orders = self.file_object.read(section,'Bonds.Orders%s'%(step_txt))
                        if isinstance(bond_orders,float) :
                                bond_orders = [bond_orders]
                        conect = {}
                        bonds = []
                        for i,(start,end) in enumerate(zip(indices[:-1],indices[1:])) :
                                if end-start > 0 :
                                        conect[i+1] = connection_table[start-1:end-1]
                                        for j in range(start-1,end-1) :
                                                bonds.append([i+1,connection_table[j],bond_orders[j]])
                        # Now correct the connection table
                        conect_sym = {}
                        for i, neighbors_i in conect.items() :
                                conect_sym[i] = neighbors_i
                                for j in neighbors_i :
                                        if not j in conect_sym.keys() :
                                                conect_sym[j] = []
                                        if j in conect.keys() : 
                                                conect_sym[j] = conect[j]
                                        if not i in conect_sym[j] :
                                                conect_sym[j].append(i)
                        conect = conect_sym     
                except KeyError :
                        pass
                return conect, bonds

        def store_mddata_for_step (self, istep) :
                """
                Store the data from the MDHistory section
                """
                if self.mddata == None : self.mddata = {}
                section = 'MDHistory'

                # First get the block info
                blocksize = self.file_object.read(section, 'blockSize')
                nblocks = self.file_object.read(section, 'nBlocks')

                # Look for the items
                sections = self.file_object.get_skeleton()
                item_keys = [kn for kn in sections[section] if 'ItemName' in kn]
                items = [self.file_object.read(section,kn) for kn in item_keys]

                # Get the data for each item
                for item in items :
                        dim = self.file_object.read(section,'%s(dim)'%(item))
                        if dim == 1 and not self.file_object.read(section,'%s(perAtom)'%(item)):
                                # Stored in block format
                                block = int(istep/blocksize)
                                pos = istep%blocksize
                                values = self.file_object.read(section,'%s(%i)'%(item,block+1))
                                if not isinstance(values,list) : values = [values]
                                self.mddata[item] = values[pos]
                        else :  
                                self.mddata[item] = self.file_object.read(section,'%s(%i)'%(item,istep+1))

        def _is_endoffile (self) :
                """
                Reads and checks If the end of file is reached.
                """
                return ('History', 'Coords(%i)'%(self.position+1)) in self.file_object

        def read_next (self, molecule=None, read=True) :
                """
                Reads the relevant info from frame self.position
                """
                if not read and not self.firsttime :
                        return self._move_cursor_without_reading()

                crd, vecs = self.read_frame (self.position,molecule)
                self.position += 1
                return crd, vecs 

        #def write_next (self, coords=None, molecule=None, cell=[0.,0.,0.], energy=0., step=None, conect=None, charges=None) :
        def write_next (self, coords=None, molecule=None, cell=[0.,0.,0.], conect=None, mddata=None) :
                """
                Write frame to next position in trajectory file

                mddata: A dictionary. It can contain the following keys:
                        'TotalEnergy', 'PotentialEnergy', 'Step', 'Velocities', 'KineticEnergy', 
                        'Charges', 'ConservedEnergy', 'Time', 'Temperature'
                """
                if isinstance(molecule,Molecule) :
                        coords, cell, elements, conect = self._read_plamsmol(molecule)
                if conect is not None :
                        if len(conect) == 0 : conect = None
                self.conect = conect
                cell = self._convert_cell(cell)

                if self.position == 0 :
                        self.write_header(coords,cell)

                step = self.position
                energy = 0.
                if mddata is not None :
                        if 'Step' in mddata :
                                step = mddata['Step']
                        if 'PotentialEnergy' in mddata :
                                energy = mddata['PotentialEnergy']

                counter = 1
                self.file_object.write('History','nEntries',self.position+1)
                self.file_object.write('History','currentEntryOpen',False)
                self._write_keydata_in_history('Step', counter, False, 1, self.position+1, step)
                counter += 1
                crd = [float(c)/bohr_to_angstrom for coord in coords for c in coord]
                self._write_keydata_in_history('Coords', counter, True, 3, self.position+1, crd)
                counter += 1
                self._write_keydata_in_history('Energy', counter, False, 1, self.position+1, energy)
                counter += 1
                if cell is not None :
                        self._write_keydata_in_history('nLatticeVectors', counter, False, 1, self.position+1, self.nvecs)
                        counter += 1
                        vecs = [float(v)/bohr_to_angstrom for vec in cell for v in vec]
                        # I should probably rethink the dimension of the lattice vectors (generalize it)
                        self._write_keydata_in_history('LatticeVectors', counter, False, [3,3], self.position+1, vecs)
                        counter += 1


                # Write the bond info
                if conect is not None :
                        # Create the index list (correct for double counting)
                        connection_table = [[at for at in conect[i+1] if at>i+1] if i+1 in conect else [] for i in range(self.ntap)]
                        numbonds = [len(neighbors) for neighbors in connection_table]
                        
                        indices = [sum(numbonds[:i])+1 for i in range(self.ntap+1)]
                        self.file_object.write('History','Bonds.Index(%i)'%(self.position+1),indices)
                        self.file_object.write('History','ItemName(%i)'%(counter),'%s'%('Bonds.Index'))
                        counter += 1

                        # Flatten the connection table
                        connection_table = [at for i in range(self.ntap) for at in connection_table[i]]
                        self.file_object.write('History','Bonds.Atoms(%i)'%(self.position+1),connection_table)
                        self.file_object.write('History','ItemName(%i)'%(counter),'%s'%('Bonds.Atoms'))
                        counter += 1

                        bond_orders = [1.0 for bond in connection_table]
                        self.file_object.write('History','Bonds.Orders(%i)'%(self.position+1),bond_orders)
                        self.file_object.write('History','ItemName(%i)'%(counter),'%s'%('Bonds.Orders'))
                        counter += 1

                if self.include_mddata and mddata is not None :
                        counter = 1
                        self.file_object.write('MDHistory','nEntries',self.position+1)
                        self.file_object.write('MDHistory','currentEntryOpen',False)
                        for key, var in mddata.items() :
                                peratom = False
                                dim = 1
                                if isinstance(var,list) :
                                        peratom = True
                                        dim = int(len(var) / self.ntap)
                                self._write_keydata_in_history(key, counter, peratom, dim, self.position+1, var, 'MDHistory')
                                counter += 1

                self.position += 1

                #if self.position%self.saving_freq == 0 : self.file_object.save()

        def _write_keydata_in_history(self, key, i, perAtom, dim, step, values, section='History') :
                """
                Write all data about a key value in KFFile
                """
                # Block code: if the data is to be written as blocks, then step and values need to be replaced.
                if section == 'MDHistory' :
                        step, values = self.get_block_info (key, perAtom, dim, step, values, section)
                                
                # The rest should be independent on format (block or individual)
                self.file_object.write(section,'%s(%i)'%(key,step),values)
                self.file_object.write(section,'ItemName(%i)'%(i),'%s'%(key))
                self.file_object.write(section,'%s(perAtom)'%(key),perAtom)
                self.file_object.write(section,'%s(dim)'%(key),dim)

        def get_block_info (self, key, perAtom, dim, step, values, section) :
                """
                If the data is to be written as blocks, then step and values need to be replaced.
                """
                if dim==1 and not perAtom :
                        try :
                                blocksize = self.file_object.read(section,'blockSize')
                        except AttributeError :
                                raise Exception ('Set include_mddata to write the MD section.')
                        iblock = int((step-1)/blocksize) + 1
                        if step%blocksize != 1 :
                                try :
                                        old_values = self.file_object.read(section,'%s(%i)'%(key,iblock))
                                        if not isinstance(old_values,list) :
                                                old_values = [old_values]
                                except AttributeError :
                                        old_values = []
                                values = old_values + [values] # Values is a scalar
                        else :
                                self.file_object.write(section,'nBlocks',iblock)
                        step = iblock
                return step, values

        def rewind (self, nframes=None) :
                """
                Rewind the file to just after the header.

                @param nframes :
                        The number of frames to rewind
                """
                self.firsttime = True
                self.position = 0

        def get_length (self) :
                """
                Get the number of steps in the file
                """
                nsteps = self.file_object.read('History', 'nEntries')
                return nsteps

        def read_last_frame (self, molecule=None) :
                """
                Reads the last geometry from the file
                """
                nsteps = self.get_length()
                crd,cell = self.read_frame(nsteps-1, molecule)
                return crd, cell
