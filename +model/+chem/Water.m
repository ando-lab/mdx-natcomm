function obj = Water()
%WATER - construct a Molecule object to represent water (hard-coded)

Name = {'O','H1','H2'}';
comp_id = {'HOH','HOH','HOH'}';
type_symbol = {'O','H','H'}';
charge = [0,0,0]';
x = [0;0.95;-.317];
y = [0;0;0];
z = [0;0;0.896];
AtomTable = table(Name,comp_id,type_symbol,charge,x,y,z);

id1 = {'O','O'}';
id2 = {'H1','H2'}';
isDouble = false(2,1);
isAromatic = false(2,1);
isCovalent = true(2,1);
BondTable = table(id1,id2,isDouble,isAromatic,isCovalent);
obj = model.chem.Molecule(AtomTable,BondTable,'name','Water','id','HOH');
end
