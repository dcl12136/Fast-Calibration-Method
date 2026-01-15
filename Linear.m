function [Obj0_T, elapsed_time] = Linear(Obj_S0, Obj_Mov, delta_L0, Dect_Std)
tic;
[M, N, K] = size(Dect_Std);
cx = floor(M/2) + 1;
cy = floor(N/2) + 1;
delta_f = 1 / (N * delta_L0);
% Compute scan positions
num_side = round(sqrt(K));
scan_positions = zeros(K, 2);
idx = 1;
for i = 1:num_side
    for j = 1:num_side
        scan_positions(idx, 1) = (i - ceil(num_side/2)) * Obj_Mov;
        scan_positions(idx, 2) = (j - ceil(num_side/2)) * Obj_Mov;
        idx = idx + 1;
    end
end
distances = sum(scan_positions.^2, 2);
[~, center_idx] = min(distances);
I_fft = zeros(M, N, K);
for k = 1:K
    I_fft(:,:,k) = fftshift(fft2(Dect_Std(:,:,k).^2));
end
I_center = I_fft(:,:,center_idx);
% Initialize spectra
O_tilde = zeros(M, N);
P_tilde = zeros(M, N);
I_dc_mean = mean(Dect_Std(:,:,center_idx).^2, 'all');
O_tilde(cx, cy) = 1.0;
P_tilde(cx, cy) = sqrt(I_dc_mean);
% Generate spiral sampling order
max_radius = min(100, floor(N/4));
spiral_order = generate_spiral_order_fast(max_radius);
num_points = size(spiral_order, 1);
% Solve spectrum point by point
solved_mask = false(M, N);
solved_mask(cx, cy) = true;
delta_x_all = scan_positions(:, 1) * delta_L0;
delta_y_all = scan_positions(:, 2) * delta_L0;
for pt_idx = 1:num_points
    di = spiral_order(pt_idx, 1);
    dj = spiral_order(pt_idx, 2);
    i_pos = cx + di;
    j_pos = cy + dj;
    i_neg = cx - di;
    j_neg = cy - dj;
    if i_pos < 1 || i_pos > M || j_pos < 1 || j_pos > N, continue; end
    if i_neg < 1 || i_neg > M || j_neg < 1 || j_neg > N, continue; end
    if solved_mask(i_pos, j_pos), continue; end
    O_ref = O_tilde(i_neg, j_neg);
    P_ref = P_tilde(i_neg, j_neg);
    if abs(O_ref) < 1e-10 || abs(P_ref) < 1e-10
        [O_ref, P_ref] = get_neighbor_ref(O_tilde, P_tilde, solved_mask, i_pos, j_pos, M, N);
    end
    zeta_center = compute_zeta_local(O_tilde, P_tilde, solved_mask, i_pos, j_pos, di, dj, M, N, 0, 0, 0);
    % Eq. 14-16
    num_frames = K - 1;
    M_mat = zeros(num_frames, 6);
    Y_vec = zeros(num_frames, 1);
    row = 0;
    for k = 1:K
        if k == center_idx, continue; end
        row = row + 1;
        delta_x = delta_x_all(k);
        delta_y = delta_y_all(k);
        zeta_k = compute_zeta_local(O_tilde, P_tilde, solved_mask, i_pos, j_pos, di, dj, M, N, delta_x, delta_y, delta_f);    
        Phi1 = 2*pi * delta_f * (-di * delta_x + dj * delta_y);
        Phi2 = 2*pi * delta_f * ( di * delta_x - dj * delta_y);
        Phi3 = 2*pi * delta_f * (2 * dj * delta_y);
        Phi4 = 2*pi * delta_f * (-(di * delta_x + dj * delta_y));  
        e1 = exp(1i * Phi1);
        e2 = exp(1i * Phi2);
        e3 = exp(1i * Phi3);
        e4 = exp(1i * Phi4);
        A1 = P_ref * (conj(zeta_k) * e1 - conj(zeta_center));
        A2 = conj(P_ref) * (zeta_k * e2 - zeta_center);
        A3 = P_ref * conj(O_ref) * (e3 - 1);
        A4 = O_ref * conj(P_ref) * (e3 - 1);
        A5 = conj(O_ref) * (zeta_k * e4 - zeta_center);
        A6 = O_ref * (conj(zeta_k) * e4 - conj(zeta_center));
        M_mat(row, :) = [A1, A2, A3, A4, A5, A6];
        Y_vec(row) = I_fft(i_pos, j_pos, k) - I_center(i_pos, j_pos);
    end
    M_mat = M_mat(1:row, :);
    Y_vec = Y_vec(1:row);
    M_real = [real(M_mat), -imag(M_mat); imag(M_mat), real(M_mat)];
    Y_real = [real(Y_vec); imag(Y_vec)];
    lambda_reg = 1e-6;
    MtM = M_real' * M_real;
    MtY = M_real' * Y_real;
    if rcond(MtM) > 1e-14
        X_real = (MtM + lambda_reg * eye(12)) \ MtY;
    else
        X_real = pinv(M_real) * Y_real;
    end
    X = X_real(1:6) + 1i * X_real(7:12);
    O_new = X(1) / max(abs(P_ref), 1e-10);
    P_new = X(5) / max(abs(conj(O_ref)), 1e-10);
    O_new = stabilize(O_new, 50);
    P_new = stabilize(P_new, 50);
    O_tilde(i_pos, j_pos) = O_new;
    P_tilde(i_pos, j_pos) = P_new;
    solved_mask(i_pos, j_pos) = true;
end
O_tilde = fill_spectrum_fast(O_tilde, solved_mask);
P_tilde = fill_spectrum_fast(P_tilde, solved_mask);
O_rec_full = ifft2(ifftshift(O_tilde));
O_amp = abs(O_rec_full);
O_phase = angle(O_rec_full);
O_phase = imgaussfilt(O_phase, 0.8);
O_amp = O_amp / max(O_amp(:));
O_amp = O_amp * 0.7 + 0.3;
O_rec_full = O_amp .* exp(1i * O_phase);
O_rec_full(~isfinite(O_rec_full)) = 1.0;
start_i = floor((M - Obj_S0) / 2) + 1;
start_j = floor((N - Obj_S0) / 2) + 1;
Obj0_T = O_rec_full(start_i:start_i+Obj_S0-1, start_j:start_j+Obj_S0-1);
elapsed_time = toc;
end 
function order = generate_spiral_order_fast(max_r)
    vals = -max_r:max_r;
    [J, I] = meshgrid(vals, vals);
    I = I(:); J = J(:);
    dists = I.^2 + J.^2;
    valid = dists > 0 & dists <= max_r^2;
    I = I(valid); J = J(valid); dists = dists(valid);
    [~, idx] = sort(dists);
    order = [I(idx), J(idx)];
end
function [O_ref, P_ref] = get_neighbor_ref(O_tilde, P_tilde, solved_mask, i, j, M, N)
    O_ref = 1.0; P_ref = 1.0;
    best = 0;
    for di = -2:2
        for dj = -2:2
            ni = i + di; nj = j + dj;
            if ni >= 1 && ni <= M && nj >= 1 && nj <= N && solved_mask(ni, nj)
                amp = abs(O_tilde(ni,nj)) + abs(P_tilde(ni,nj));
                if amp > best
                    best = amp;
                    O_ref = O_tilde(ni,nj);
                    P_ref = P_tilde(ni,nj);
                end
            end
        end
    end
end
function zeta = compute_zeta_local(O_tilde, P_tilde, solved_mask, i, j, di, dj, M, N, dx, dy, df)
    zeta = 0;
    for delta_i = -2:2
        for delta_j = -2:2
            if delta_i == 0 && delta_j == 0, continue; end
            ni = i + delta_i;
            nj = j + delta_j;
            if ni < 1 || ni > M || nj < 1 || nj > N, continue; end
            if ~solved_mask(ni, nj), continue; end
            neighbor_di = di + delta_i;
            neighbor_dj = dj + delta_j;
            if (neighbor_di == di && neighbor_dj == dj) || (neighbor_di == -di && neighbor_dj == -dj)
                continue;
            end
            O_val = O_tilde(ni, nj);
            P_val = P_tilde(ni, nj);
            if dx == 0 && dy == 0 && df == 0
                zeta = zeta + O_val * P_val;
            else
                phase = 2 * pi * df * (neighbor_di * dx + neighbor_dj * dy);
                zeta = zeta + O_val * P_val * exp(1i * phase);
            end
        end
    end
end
function val = stabilize(val, max_amp)
    if ~isfinite(val)
        val = 0;
    elseif abs(val) > max_amp
        val = val / abs(val) * max_amp;
    end
end
function spectrum = fill_spectrum_fast(spectrum, solved_mask)
    kernel = ones(3,3)/8; kernel(2,2) = 0;
    for iter = 1:3
        filled = conv2(spectrum, kernel, 'same');
        need_fill = ~solved_mask & abs(spectrum) < 1e-10;
        spectrum(need_fill) = filled(need_fill) * 0.8;
    end
end