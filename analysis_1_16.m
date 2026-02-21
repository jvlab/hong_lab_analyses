X = tdata{:,:};
X = X-repmat(mean(X,2),[1,size(X,2)]);

X0 = X(1:8,:);
X1 = X(9:23,:);
numGlom = size(X,2);
% remove the column means
%X0_rm = X0-repmat(mean(X0),[8 1]);
%X1_rm = X1-repmat(mean(x1),[15 1]);


X0_rm = X0;
X1_rm = X1;
X0norm = diag(X0_rm*X0_rm');
X0norm = 1./X0norm;
X0norm = sqrt(X0norm);
X0norm = diag(X0norm);
X0normalized = X0norm * X0;
X1norm = diag(X1_rm*X1_rm');
X1norm = 1./X1norm;
X1norm = sqrt(X1norm);
X1norm = diag(X1norm);
X1normalized = X1norm * X1;
    

[U,S,V] = svd(X0,'econ'); % matlab passes v, not vt like everything else.

A0 = X0*V;
A1 = X0normalized*V;
B0 = X1*V;
B1 = X1normalized*V;

figure;
% ms


msm = [A1(1:2,:);B1(1:4,:)]';
%figure;
subplot(2,3,1)
msb = bar(msm);
msb(1).FaceColor=[0 0 1];
msb(2).FaceColor=[0 0 1];
for it=3:6
    msb(it).FaceColor=[1 0 0];
end
title('ms(-3.0)')

vam = [A1(3:4,:);B1(5:6,:)]';
%figure;
subplot(2,3,2)
vab = bar(vam);
for it=1:2
    vab(it).FaceColor=[0 0 1];
    vab(it+2).FaceColor=[1 0 0];
end
title('va(-3.0)');

%2h
twohm = [A1(5,:);B1(7:10,:)]';
subplot(2,3,3)
twohmb = bar(twohm);
twohmb(1).FaceColor=[0 0 1];
for it=2:5
    twohmb(it).FaceColor = [1 0 0];
end
title('2h(-6.0)');

%p-cre
pcrem = [A1(7,:);B1(11:13,:)]';
subplot(2,3,4)
pcreb = bar(pcrem);
pcreb(1).FaceColor=[0 0 1];
for it=2:4
    pcreb(it).FaceColor=[1 0 0];
end
title('p-cre');

%ACV

ACVm = [A1(7,:);B1(14,:)]';
subplot(2,3,5)
ACVb = bar(ACVm);
ACVb(1).FaceColor = [0 0 1];
ACVb(2).FaceColor = [1 0 0];
title('ACV')

%Iaa

Iaam = [A1(8,:);B1(15,:)]';
subplot(2,3,6)
Iaab = bar(Iaam);
Iaab(1).FaceColor=[0 0 1];
Iaab(2).FaceColor=[1 0 0];
title('Iaa')