X0 = standard{:,:};
X0 = X0-repmat(mean(X0,2),[1 size(X0,2)]);
X1 = diag_repeat{:,:};
X1 = X1-repmat(mean(X1,2),[1 size(X1,2)]);

numGlom = size(X0,2);

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
indx = contains(standard.Properties.RowNames,'ms(-3.0)');
msm = [A1(indx,:);B1(1:4,:)]';
%msm = [A1(1:2,:);B1(1:4,:)]';
%figure;
subplot(2,3,1)
msb = bar(msm);
msb(1).FaceColor=[0 0 1];
msb(2).FaceColor=[0 0 1];
for it=3:6
    msb(it).FaceColor=[1 0 0];
end
title('ms(-3.0)')
%
indx = contains(standard.Properties.RowNames,'va(-3.0)');
vam = [A1(indx,:);B1(5:6,:)]';
%figure;
subplot(2,3,2)
vab = bar(vam);
for it=1:2
    vab(it).FaceColor=[0 0 1];
    vab(it+2).FaceColor=[1 0 0];
end
title('va(-3.0)');
%
%2h
indx = startsWith(standard.Properties.RowNames,'2h(-6.0)');
twohm = [A1(indx,:);B1(7:10,:)]';
subplot(2,3,3)
twohmb = bar(twohm);
twohmb(1).FaceColor=[0 0 1];
for it=2:5
    twohmb(it).FaceColor = [1 0 0];
end
title('2h(-6.0)');

% 
%p-cre
indx = startsWith(standard.Properties.RowNames,'p-cre(-3.0)');
pcrem = [A1(indx,:);B1(11:14,:)]';
subplot(2,3,4)
pcreb = bar(pcrem);
pcreb(1).FaceColor=[0 0 1];
for it=2:5
    pcreb(it).FaceColor=[1 0 0];
end
title('p-cre');

%ACV
%
indx = startsWith(standard.Properties.RowNames,'ACV(0.0)');

ACVm = [A1(indx,:);B1(15,:)]';
subplot(2,3,5)
ACVb = bar(ACVm);
ACVb(1).FaceColor = [0 0 1];
ACVb(2).FaceColor = [1 0 0];
title('ACV')

%Iaa
indx = startsWith(standard.Properties.RowNames,'IaA(-1.0)');

Iaam = [A1(indx,:);B1(16,:)]';
subplot(2,3,6)
Iaab = bar(Iaam);
Iaab(1).FaceColor=[0 0 1];
Iaab(2).FaceColor=[1 0 0];
title('Iaa')
%}