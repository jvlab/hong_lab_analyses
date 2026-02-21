% Check the relationship between normalized odorants at different
% concentrations.
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
%1p3ol
figure;

odors = startsWith(standard.Properties.RowNames,'1p3ol');

A_1p3ol = A1(odors,:)';

subplot(2,3,1);
msm=bar(A_1p3ol);
title('1p3ol');
% 2-but
odors = startsWith(standard.Properties.RowNames,'2-but(');

A_2but = A1(odors,:)';
subplot(2,3,2);
ms2=bar(A_2but);
title('2-but');

% 2h
odors = startsWith(standard.Properties.RowNames,'2h(');

A_2h = A1(odors,:)';
subplot(2,3,3);
ms3=bar(A_2h);
title('2h');

% p-cre
odors = startsWith(standard.Properties.RowNames,'p-cre(');

A_pcre = A1(odors,:)';
subplot(2,3,4);
ms4=bar(A_pcre);
title('p-cre');

% banana
odors = startsWith(standard.Properties.RowNames,'ban(');

A_ban = A1(odors,:)';
subplot(2,3,5);
ms5=bar(A_2but);
title('banana');

% t2h
odors = startsWith(standard.Properties.RowNames,'peara(');

A_peara = A1(odors,:)';
subplot(2,3,6);
ms6=bar(A_peara);
title('peara')