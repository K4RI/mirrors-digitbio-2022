function cca_obj_test=LCCA(file, indX, indY)

warning('off','all');

data = dlmread(file);

pkg load statistics

x=data(:,indX);
y=data(:,indY);

total_number_of_samples=size(x,1);
permutation=randperm(size(x,1))';         %reandomly permute the observations
sampled_samples=total_number_of_samples;  %how many observations to sample from the data

%train_size=floor(0.7*sampled_samples);
%ev_size=floor(0.15*sampled_samples);
%test_size=floor(0.15*sampled_samples);

%x_train=x(permutation(1:train_size),:);
%y_train=y(permutation(1:train_size),:);

%x_ev=x(permutation(train_size+1:train_size+ev_size),:);
%y_ev=y(permutation(train_size+1:train_size+ev_size),:);
%x_test=x(permutation(train_size+ev_size+1:sampled_samples),:);
%y_test=y(permutation(train_size+ev_size+1:sampled_samples),:);

%x_train_ev=[x_train; x_ev];
%y_train_ev=[y_train; y_ev];

x_train = x;
y_train = y;
x_ev = x;
y_ev = y;
x_test = x;
y_test = y;
x_train_ev = x;
y_train_ev = y;

[a,b,r] = canoncorr(x_train_ev,y_train_ev);
cca_obj_train=sum(r);
x_cca_test=x_test*a;
y_cca_test=y_test*b;
cov_mat=cov([x_cca_test,y_cca_test]);
cov_mat=cov_mat(1:size(x_cca_test,2),size(x_cca_test,2)+1:end);
r_mat=(cov_mat.^2)./(cov(x_cca_test).*cov(y_cca_test));
cca_obj_test=sum(sqrt(diag(r_mat)));

end