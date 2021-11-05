function [cams_cluster] = create_graph_patches(A,opts)
tic;
[I,J]=find(tril(A,-1));
ncams=size(A,1);
B = adj2inc(A);
initial_datum = opts.initial_datum;
min_legs = opts.min_legs;
max_size= floor(opts.max_size_coefficient*3.7*sqrt(ncams));
min_size= floor(max_size/2);
% Clustering
ngroups=ceil(ncams/min_size);
groups = spectral_shi_malik(A,ngroups,initial_datum);
cams_cluster=cell(ngroups,1);
for k=1:ngroups
    cams_cluster{k} = find(groups==k);
end
t0 = toc;
end