function [B0_motion] = FindCharacteristicMotionMap(B0_SVD_struct,coef_bias,coef_dc,coef_R2,coef_T_rigid,num_vectors,n_ix,p_ix)

T_rigid_motion_p = B0_SVD_struct.T_rigid_p;
T_rigid_motion_n = B0_SVD_struct.T_rigid_n;
weight_motionp = B0_SVD_struct.weight_p;
weight_motionn = B0_SVD_struct.weight_n;
B0_SVD_motionp = B0_SVD_struct.B0_SVD_p;
B0_SVD_motionn = B0_SVD_struct.B0_SVD_n;
coefs_motionp = B0_SVD_struct.coefs_p;
coefs_motionn = B0_SVD_struct.coefs_n;

figure;
row = 1;
for i = 1:num_vectors
    subplot(num_vectors,2,row);
    plot(T_rigid_motion_n(coef_T_rigid,:),weight_motionn(:,i),'o','linewidth',2); xlabel('FH translation (mm)'); ylabel('Vector weight (-)'); hold on;
    x_axis_motion = linspace(min(T_rigid_motion_n(coef_T_rigid,:)),max(T_rigid_motion_n(coef_T_rigid,:)),100);
    plot(x_axis_motion,coefs_motionn(coef_bias,i)+x_axis_motion*coefs_motionn(coef_dc,i)); title(['(',num2str(i),')']); hold off;
    subplot(num_vectors,2,row+1);
    plot(T_rigid_motion_p(coef_T_rigid,:),weight_motionp(:,i),'o','linewidth',2); xlabel('FH translation (mm)'); ylabel('Vector weight (-)'); hold on;
    x_axis_motion = linspace(min(T_rigid_motion_p(coef_T_rigid,:)),max(T_rigid_motion_p(coef_T_rigid,:)),100);
    plot(x_axis_motion,coefs_motionp(coef_bias,i)+x_axis_motion*coefs_motionp(coef_dc,i)); title(['(',num2str(i),')']); hold off;
    row = row + 2;
end

% For basis vector, the vector with highecst R2 is chosen
if n_ix == 0
    [~,motionn_ix] = max(coefs_motionn(coef_R2,:)); % R2 w.r.t. negative FH translations
else
    motionn_ix = n_ix;
end
if p_ix == 0
    [~,motionp_ix] = max(coefs_motionp(coef_R2,:)); % R2 w.r.t. positive FH translations
else
    motionp_ix = p_ix;
end
disp(motionn_ix);

% Make sure the B0 map entails a 1 degree rotation in exclusively 1 plane
B0_motion = zeros([size(B0_SVD_motionn,1),size(B0_SVD_motionn,2),2]);
B0_motion(:,:,1) = 1.00*B0_SVD_motionn(:,:,motionn_ix)*(coefs_motionn(coef_bias,motionn_ix)+1.0*coefs_motionn(coef_dc,motionn_ix)); 
B0_motion(:,:,2) = 1.00*B0_SVD_motionp(:,:,motionp_ix)*(coefs_motionp(coef_bias,motionp_ix)+1.0*coefs_motionp(coef_dc,motionp_ix)); 

figure;
subplot(1,2,1);
imagesc(2*pi*B0_motion(:,:,1), [-50 50]); title('Negative FH translation'); colorbar; axis image; axis off;
subplot(1,2,2);
imagesc(2*pi*B0_motion(:,:,2), [-50 50]); title('Positive FH translation'); colorbar; axis image; axis off;

end