function M = convert2Molecule(mmcif)

assert(numel(mmcif)==1,'only one molecule allowed at a time');
x = str2double(mmcif.chem_comp_atom.model_Cartn_x);
y = str2double(mmcif.chem_comp_atom.model_Cartn_y);
z = str2double(mmcif.chem_comp_atom.model_Cartn_z);
comp_id = strrep(mmcif.chem_comp_atom.comp_id,'"','');
charge = str2double(mmcif.chem_comp_atom.charge);
type_symbol = strrep(mmcif.chem_comp_atom.type_symbol,'"','');

Names = strrep(mmcif.chem_comp_atom.atom_id,'"','');

if isfield(mmcif,'chem_comp_bond')

    isDouble = strcmpi(mmcif.chem_comp_bond.value_order,'doub');
    isAromatic = strcmpi(mmcif.chem_comp_bond.pdbx_aromatic_flag,'y');
    isCovalent = true(size(isDouble));
    
    atom_id_1 = strrep(mmcif.chem_comp_bond.atom_id_1,'"','');
    atom_id_2 = strrep(mmcif.chem_comp_bond.atom_id_2,'"','');
    BondTable = table(atom_id_1,atom_id_2,isDouble,isAromatic,isCovalent);
else

    BondTable = table([],[],[],[],[],'variableNames',...
        {'atom_id_1','atom_id_2','isDouble','isAromatic','isCovalent'});
end
AtomTable = table(Names,comp_id,type_symbol,charge,x,y,z);

mol_name = strrep(mmcif.chem_comp.name,'"','');
mol_id = strrep(mmcif.chem_comp.id,'"','');

M = model.chem.Molecule(AtomTable,BondTable,'name',mol_name,'id',mol_id);

end