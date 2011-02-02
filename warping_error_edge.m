function [warping_error, warped_labels, nonsimple_classify]=warping_error_edge(label, proposal, mask, binary_threshold, radius, fg_conn)

if(~exist('mask') | isempty(mask))
	mask=ones(size(label));
end

if(~exist('binary_threshold') | isempty(binary_threshold))
	binary_threshold=0.5;
end

if(~exist('fg_conn') | isempty(binary_threshold))
	fg_conn=6;
end

if(~exist('radius') | isempty(radius))
	radius=Inf;
end

if radius==Inf
	radius=0;
end

%fprintf('Stage 1\n');
warped_labels=simple_edge_warp_3d_mex(uint32(label), single(proposal), logical(mask~=0), single(binary_threshold),uint32(fg_conn));
%fprintf('Stage 2\n');
error_mask=((proposal>binary_threshold)~=(warped_labels~=0)).*max(mask,[],5);
%fprintf('Stage 3\n');
nonsimple_classify=non_simple_3d_edge_classify_mex_parallel(uint32(warped_labels), logical(error_mask~=0), uint32(radius), uint32(fg_conn));
%fprintf('Stage 4\n');
warping_error=nnz(error_mask)/nnz(mask);
%fprintf('Stage 5\n');
