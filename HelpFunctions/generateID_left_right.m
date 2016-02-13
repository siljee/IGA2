function ID = generateID_left_right(np_xi, np_eta, isLeftDirichlet, isRightDirichlet)

Dbnd = 0;                   % Number of dirichlet bnd (left, right or both)

if ~isempty(isLeftDirichlet)
   Dbnd = Dbnd +1;
end
if ~isempty(isRightDirichlet)
   Dbnd = Dbnd +1;
end


ID = [];
for i = 1:np_eta
    ID = [ID,    isLeftDirichlet,(i-1)*(np_xi-Dbnd)+1:i*(np_xi-Dbnd),isRightDirichlet];
end