function IEN = generate_IEN(knotVec_xi, knotVec_eta, p_xi, p_eta)

np_xi = length(knotVec_xi) - p_xi -1;                   % Number of control points in xi direction
np_eta = length(knotVec_eta) - p_eta -1;                % Number of control points in eta direction

% To account for multplicazity in knot vectors
mult_xi = abs((knotVec_xi(p_xi+1:end-p_xi-1)-knotVec_xi(p_xi+2:end-p_xi)~= 0) - 1);
mult_eta = abs((knotVec_eta(p_eta+1:end-p_eta-1)-knotVec_eta(p_eta+2:end-p_eta)~= 0) - 1);

el_xi = np_xi-p_xi-sum(mult_xi);
el_eta = np_eta-p_eta-sum(mult_eta);


IEN=zeros(np_xi*np_eta,(p_xi+1)*(p_eta+1));

row = [];
col = [];
col1 = [];

round = 0;
for i = 1:p_eta+1
    row = [row, (1:p_xi+1) + round*(np_xi)];
    round = round + 1;
end

round = 0;
for i = mult_xi%np_eta
    if i == 0
        col1 = [col1, round];%[col, (0:np_xi-1-p_xi) + round*np_xi]
    end
    round = round + 1;
end

round = 0;
for i = mult_eta
    if i == 0
        col = [col col1+round*np_xi];
    end
    round = round +1;
end

row = ones(el_xi*el_eta,1)*row;
col = col'*ones(1,(p_xi+1)*(p_eta+1));
IEN = row + col;
