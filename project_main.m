%%
clear;
%% Section 1: Functional Verification of Your MATLAB Reference Model
load CEN571_class_project_data.mat;

%%

num_matrix = size(A,3);
m = size(A,1);
n = size(A,2);
total_operations = num_matrix*2*m*n^2;

for i = 1:num_matrix
   [Q(:,:,i), R(:,:,i)] = myQR (A(:,:,i));
end

MSE_ref = getMSE(Q, R, A)

if MSE_ref < 1e-28
    disp('Functional Verification is Successful. ^_^');
else
    disp('Functional Verification is Failed. x_x');
end

%% Section 2: Caculate and Set Your Simulation Stop Time
% Calculate the simulation stop time as the time when the last Q and R
% results are completely produced. It will be used to calculate the
% computational throughput of your accelerator design. So make sure to use
% the minimum time that can produce all the Q and R results. 

input_sample_freq = 1e6; % Input Sampling Frequency is 1 MHz
sample_time = 1 / input_sample_freq; % You must set the sample time of your accelerator subsystem as the sample_time virable.
% iteration_bound = 5;
% Note : The feedforward latency described the variable
% 'feedforward_latency' is the total latency of the fixed point design and
% its used to calculate the simulation stop time for the same
feedforward_latency = 37;%20 % 32
interleaving_factor = 1;


num_samples_test = num_matrix;

time = [0:num_samples_test-1]'*sample_time;

simulation_stop_time = sample_time * (num_samples_test + feedforward_latency/interleaving_factor); % in second, caculate and set the stop time of your simulation as the simulation_stop_time variable.
cycle_runtime = simulation_stop_time / sample_time; %The cycle runtime of your design will be calculated as simulation_stop_time / sample_time.

%% Section 3: Simulate Your Double-Precision Accelerator Subsystem and Collect Results
% Design your own accelerator subsystem and testbench. Simulate your 
% double-precision accelerator subsystem with the testing data A and the 
% sample_time and simulation_stop_time defined above. Collect the
% output results and pack your results into Q_hdl_flp and R_hdl_flp, each 
% being a 5x5x1000 tensor. 

Array_A_Input = reshape(A,[m*n*num_matrix,1]);

Array_A_Input = A(:,:,1:num_samples_test);

% At.time = time;
% At.signals.values = Array_A_Input ;
% At.signals.dimensions = [5 5];


for i = 1:5
   AtI(i).time = time;
   AtI(i).signals.values = Array_A_Input(:,i,:);
   AtI(i).signals.dimensions = [5 1];
end




%% Section 3.1:  Now computing the outputs  indexes and concatenation indexes (Pre Processing step for indexing the output values)

%Note : The variable 'feedforward_latency_flp' is the latency within the
%Floating point design and this value is used to index the output values
%obtained after running the Floating point design.
feedforward_latency_flp = 20;

system_latency = feedforward_latency_flp/5;

index1 = [];
index2 = [];
index3 = [];
index4 = [];
index5 = [];

for i = 1:num_samples_test
    
index1 = [index1 (system_latency+i)];
index2 = [index2 (system_latency*2 + i)];
index3 = [index3 (system_latency*3 + i)];
index4 = [index4 (system_latency*4 + i)];
index5 = [index5 (feedforward_latency_flp + i)];


end

%% Section 3.2: Collecting floating point results 
Q_hdl_flp = zeros(5,5,num_samples_test);
R_hdl_flp = zeros(5,5,num_samples_test);

for i=1:num_samples_test
   
    q1 = out.q1.signals.values(:,:,index1(i));
    q2 = out.q2.signals.values(:,:,index2(i));
    q3 = out.q3.signals.values(:,:,index3(i));
    q4 = out.q4.signals.values(:,:,index4(i));
    q5 = out.q5.signals.values(:,:,index5(i));
    Q = [q1,q2,q3,q4,q5];
    
    Q_hdl_flp(:,:,i) = Q;
    
    
    r1 = out.r1.signals.values(:,:,index1(i));
    r2 = out.r2.signals.values(:,:,index2(i));
    r3 = out.r3.signals.values(:,:,index3(i));
    r4 = out.r4.signals.values(:,:,index4(i));
    r5 = out.r5.signals.values(:,:,index5(i));
    R = [r1;r2;r3;r4;r5];
    
    R_hdl_flp(:,:,i) = R;
    
end




%% Section 4: Functional Verification of Your Double-Precision Accelerator Subsystem

MSE_flp = getMSE(Q_hdl_flp, R_hdl_flp, Array_A_Input)

if MSE_flp < 1e-28
    disp('Functional Verification is Successful. ^_^');
else
    disp('Functional Verification is Failed. x_x');
end

%% Section 5: Data Type Optimization
% Optimize the data types of your design with a relative tolerence of 1e-4
% applied to each output of your design. 

model_name = 'QRF_fp';
system_name = 'QR_factorization_engine';
QR_Engine = [model_name '/' system_name];
open_system(model_name);

opt = fxpOptimizationOptions('AllowableWordLengths', 8:72, 'UseParallel', false);

% define the stopping criterion in terms of the maximum number of
% iterations
opt.MaxIterations = 100;

% define an absolute tolerance at the output of the system to constrain 
% the differences between the original output values of the system 
% and that based on the optimized fixed-point data types.
tol1 = 1e-1;
tol2 = 1e-1;

addTolerance(opt, [model_name, '/Q_Out1'], 1, 'AbsTol', tol1);
addTolerance(opt, [model_name, '/Q_Out2'], 1, 'AbsTol', tol1);
addTolerance(opt, [model_name, '/Q_Out3'], 1, 'AbsTol', tol1);
addTolerance(opt, [model_name, '/Q_Out4'], 1, 'AbsTol', tol1);
addTolerance(opt, [model_name, '/Q_Out5'], 1, 'AbsTol', tol1);


addTolerance(opt, [model_name, '/R_Out1'], 1, 'AbsTol', tol2);
addTolerance(opt, [model_name, '/R_Out2'], 1, 'AbsTol', tol2);
addTolerance(opt, [model_name, '/R_Out3'], 1, 'AbsTol', tol2);
addTolerance(opt, [model_name, '/R_Out4'], 1, 'AbsTol', tol2);
addTolerance(opt, [model_name, '/R_Out5'], 1, 'AbsTol', tol2);


% Use fxpopt function to run the optimization
result = fxpopt(model_name, QR_Engine, opt);
% Launch Simulation Data Inspector
explore(result);


%% Section 6.0:  Now computing the output indexes and concatenation Index for the fixed point version (Pre processing for output collection of Fixed point design)
%Note : The fixed point design is pipelined in a certain way 
% Hence we have 3 variables that represent that pipelining 
% ffwd_q is the delay for producing the first output column of Q for each
% 5X5 matrix as input 
% ffwd_r is the delay for producing the first output row of R for each 5X5
% matrix as input 
% cumulative_ff is the delay between two Q and R outputs, i.e., two columns
% of Q or two rows of R at the output 
% These values are used to compute proper indexes of outputs so that we can
% collect the results properly 

ffwd_q = 4;
ffwd_r = 6;
cumulative_ff = 8;

index1_q = [];
index2_q = [];
index3_q = [];
index4_q = [];
index5_q = [];

index1_r = [];
index2_r = [];
index3_r = [];
index4_r = [];
index5_r = [];


for i = 1:num_samples_test
    
index1_q = [index1_q (ffwd_q+i)];
index1_r = [index1_r (ffwd_r+i)];

index2_q = [index2_q (ffwd_q + cumulative_ff*1 + i)];
index2_r = [index2_r (ffwd_r + cumulative_ff*1 + i)];

index3_q = [index3_q (ffwd_q + cumulative_ff*2 + i)];
index3_r = [index3_r (ffwd_r + cumulative_ff*2 + i)];

index4_q = [index4_q (ffwd_q + cumulative_ff*3 + i)];
index4_r = [index4_r (ffwd_r + cumulative_ff*3 + i)];

index5_q = [index5_q (ffwd_q + cumulative_ff*4 + i)];
index5_r = [index5_r (ffwd_r + cumulative_ff*4 + i)];


end

%% Section 6: Simulate Your Data-Type-Optimized Accelerator Subsystem and Collect Results
% Simulate the data-type-optimized accelerator subsystem with the testing 
% data A and the sample_time and simulation_stop_time defined above. 
% Collect the output results and pack your results into Q_hdl_fp and 
% R_hdl_fp, each being an 5x5x1000 tensor. Calculate the new MSE at the 
% reduced numerical precision. 
Q_hdl_fp = zeros(5,5,num_samples_test);
R_hdl_fp = zeros(5,5,num_samples_test);

for i=1:num_samples_test
   
    q1 = out.q1.signals.values(:,:,index1_q(i));
    q2 = out.q2.signals.values(:,:,index2_q(i));
    q3 = out.q3.signals.values(:,:,index3_q(i));
    q4 = out.q4.signals.values(:,:,index4_q(i));
    q5 = out.q5.signals.values(:,:,index5_q(i));
    Q = [q1,q2,q3,q4,q5];
    
    Q_hdl_fp(:,:,i) = Q;
    
    
    r1 = out.r1.signals.values(:,:,index1_r(i));
    r2 = out.r2.signals.values(:,:,index2_r(i));
    r3 = out.r3.signals.values(:,:,index3_r(i));
    r4 = out.r4.signals.values(:,:,index4_r(i));
    r5 = out.r5.signals.values(:,:,index5_r(i));
    R = [r1;r2;r3;r4;r5];
    
    R_hdl_fp(:,:,i) = R;
    
end


MSE_fp = getMSE(Q_hdl_fp, R_hdl_fp, Array_A_Input)

%% Section 7: FPGA Mapping
% Map the data-type-optimized design into a Xilinx XC7A100T-CSG324C (Speed -1)
% FPGA target. Check the utilization report details and report the FPGA 
% resource utilizations. Make sure the utilization of each type of FPGA 
% resources is <=100%. Then, check the timing report details to calculate 
% and report the maximum clock frequency your accelerator design can run at.
% Make sure you have a timing closure. 

% Report the utilization percentage, e.g. 0.01 means 1%. Must be <= 1 (100%).
slice_LUTs_util = 0.4862;
slice_Registers_util = 0.0895;
DSPs_util = 1.00; 
Block_RAM_Tile_util = 0.1522;

% Report the target frequency and the worst negative slack (WNS)
target_frequency = 2; % MHz
worst_negative_slack = 287.830 % ns, WNS reported by synthesis. Must be >= 0. 


actual_clock_period = 1e3 / target_frequency - worst_negative_slack; %ns
actual_max_frequency = 1e3 / actual_clock_period %MHz

%% Section 8: Calculate Computational Throughput
% The computational througput of your design
% will be calculated as total_operations / (cycle_runtime / actual_max_frequency ). 

actual_runtime = cycle_runtime / (actual_max_frequency*1e6)
computational_throughput = total_operations / actual_runtime / 1e9 % GOPS

