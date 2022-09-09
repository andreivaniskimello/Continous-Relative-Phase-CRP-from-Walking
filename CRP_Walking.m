%Walking CRP from Segments
%by André Ivaniski Mello (andreivaniskimello@gmail.com)
%Reference: Lamb & Stockl, 2014, doi: 10.1016/j.clinbiomech.2014.03.008

%INPUT: 
% Angles curves from body segments (0-100% from stride cycle)
% Subject Name and Condition

% INSERT ALSO:

% Butterworth Filter Parameters desired


%OUTPUT:

%Excel file with: Segments pair CRP



    % CRP Calculation Steps (Lamb & Stockl, 2014)

        % 1. centering the amplitude of the data around zero (Eq. (11))
        % 2. transforming each signal into an analytic signal using the Hilbert transform (Eq. (7))
        % 3. calculating the phase angles for each signal (Eq. (8))
        % 4. calculating the continuous relative phase (Eq. (9))

    % Step 1)
    % x_centered = x - min(x) - (max(x) - min(x))/2;

    % Step 2)
    % H = H(x);
    % x = hilbert(xr)
    % imx = imag(x) %The imaginary part of x is the Hilbert transform of xr, and the real part is xr itself.
        %To CRP, use only the imaginary part of x

    % Step 3)
    % Phase_Angle = arctan(H(x) / x);

    % Step 4)
    % CRP = Phase_Angle_Segment_1_Proximal - Phase_Angle_Segment_2_Distal;


clear
close all
clc




%% SECTION 01
% INFORMATION FOR ROUTINE %

%Insert Input Files Folder Path here: this is the folder where all Kinematic individual
%files are contained
Input_File_Folder = 'C:\Users\andre\Documents\Andre\Pesquisa\Artigos para Publicar\Mestrado - Biomechanics Shallow Water Walking\Data\Data In - Angular Curves\';

%Insert Output File Path here
Output_File_Path = 'C:\Users\andre\Documents\Andre\Pesquisa\Artigos para Publicar\Mestrado - Biomechanics Shallow Water Walking\Data\Data Out\'; 

% Information File Data Path here (is only ONE file for all Kinematic individual data files)
Input_Files_Information = "C:\Users\andre\Documents\Andre\Pesquisa\Artigos para Publicar\Mestrado - Biomechanics Shallow Water Walking\Data\Angular_Curves_Information.xlsx";

Input_Matrix_FilesNames = readtable(Input_Files_Information);
Input_Matrix_FilesNames = table2array(Input_Matrix_FilesNames(:,1));

% Here the numerical data from the from TDTO File is imported
Input_File_Data = readmatrix(Input_Files_Information);



%Depth code
    % 1 = Knee
    % 2 = Hip
    % 3 = Umbilical/Xiphoid
    
%%

    
for k =  1:length(Input_File_Data)
    
    % Load kinematic data

        % Here the File Name from the Path selected is read.
    
    File_Name_Individual = char(Input_Matrix_FilesNames(k));
    File_Name_Path_Full = [Input_File_Folder File_Name_Individual];
    fprintf(1, 'Now reading %s\n', File_Name_Individual);
    
    % LOAD Input Data (One file per iteration)
    Delimiter_Columns = '\t';
    Header_Lines = 5;
    Input_Data = importdata(File_Name_Path_Full, Delimiter_Columns, Header_Lines);
    Data = Input_Data.data;
    
    %Extract the depth and speed condition
    Depth = Input_File_Data(k,2);
    Speed = Input_File_Data(k,3);
    
    if Depth == 1 
    
        % Knee Depth
        
        % Extract the Angles curves from each segment and joint. 5 columns correspond to the 5 strides.
        
        Foot_Angle = Data(:,1:5);
        Shank_Angle = Data(:,7:11);

        Ankle_Joint_Angle = Data(:,13:17);

                %Separating per Strides
        Foot_Angle_Stride_1 = Foot_Angle(:,1);
        Foot_Angle_Stride_2 = Foot_Angle(:,2);
        Foot_Angle_Stride_3 = Foot_Angle(:,3);
        Foot_Angle_Stride_4 = Foot_Angle(:,4);
        Foot_Angle_Stride_5 = Foot_Angle(:,5);

        Shank_Angle_Stride_1 = Shank_Angle(:,1);
        Shank_Angle_Stride_2 = Shank_Angle(:,2);
        Shank_Angle_Stride_3 = Shank_Angle(:,3);
        Shank_Angle_Stride_4 = Shank_Angle(:,4);
        Shank_Angle_Stride_5 = Shank_Angle(:,5);

        Ankle_Joint_Angle_Stride_1 = Ankle_Joint_Angle(:,1);
        Ankle_Joint_Angle_Stride_2 = Ankle_Joint_Angle(:,2);
        Ankle_Joint_Angle_Stride_3 = Ankle_Joint_Angle(:,3);
        Ankle_Joint_Angle_Stride_4 = Ankle_Joint_Angle(:,4);
        Ankle_Joint_Angle_Stride_5 = Ankle_Joint_Angle(:,5);

        % CRP Calculation

            % Step 1 Centering angle
            for i = 1:length(Foot_Angle_Stride_1)
                Foot_Angle_Stride_1_Centered(i) =  Foot_Angle_Stride_1(i) - min(Foot_Angle_Stride_1) -  (( max(Foot_Angle_Stride_1) - min(Foot_Angle_Stride_1) )  /2);
                Foot_Angle_Stride_2_Centered(i) =  Foot_Angle_Stride_2(i) - min(Foot_Angle_Stride_2) -  (( max(Foot_Angle_Stride_2) - min(Foot_Angle_Stride_2) )  /2);
                Foot_Angle_Stride_3_Centered(i) =  Foot_Angle_Stride_3(i) - min(Foot_Angle_Stride_3) -  (( max(Foot_Angle_Stride_3) - min(Foot_Angle_Stride_3) )  /2);
                Foot_Angle_Stride_4_Centered(i) =  Foot_Angle_Stride_4(i) - min(Foot_Angle_Stride_4) -  (( max(Foot_Angle_Stride_4) - min(Foot_Angle_Stride_4) )  /2);
                Foot_Angle_Stride_5_Centered(i) =  Foot_Angle_Stride_5(i) - min(Foot_Angle_Stride_5) -  (( max(Foot_Angle_Stride_5) - min(Foot_Angle_Stride_5) )  /2);
                
                Shank_Angle_Stride_1_Centered(i) =  Shank_Angle_Stride_1(i) - min(Shank_Angle_Stride_1) -  (( max(Shank_Angle_Stride_1) - min(Shank_Angle_Stride_1) )  /2);
                Shank_Angle_Stride_2_Centered(i) =  Shank_Angle_Stride_2(i) - min(Shank_Angle_Stride_2) -  (( max(Shank_Angle_Stride_2) - min(Shank_Angle_Stride_2) )  /2);
                Shank_Angle_Stride_3_Centered(i) =  Shank_Angle_Stride_3(i) - min(Shank_Angle_Stride_3) -  (( max(Shank_Angle_Stride_3) - min(Shank_Angle_Stride_3) )  /2);
                Shank_Angle_Stride_4_Centered(i) =  Shank_Angle_Stride_4(i) - min(Shank_Angle_Stride_4) -  (( max(Shank_Angle_Stride_4) - min(Shank_Angle_Stride_4) )  /2);
                Shank_Angle_Stride_5_Centered(i) =  Shank_Angle_Stride_5(i) - min(Shank_Angle_Stride_5) -  (( max(Shank_Angle_Stride_5) - min(Shank_Angle_Stride_5) )  /2);
                
            end

            % Step 2 Hilbert transform
            Foot_Angle_Stride_1_Hilbert = hilbert(Foot_Angle_Stride_1_Centered);
            Foot_Angle_Stride_2_Hilbert = hilbert(Foot_Angle_Stride_2_Centered);
            Foot_Angle_Stride_3_Hilbert = hilbert(Foot_Angle_Stride_3_Centered);
            Foot_Angle_Stride_4_Hilbert = hilbert(Foot_Angle_Stride_4_Centered);
            Foot_Angle_Stride_5_Hilbert = hilbert(Foot_Angle_Stride_5_Centered);
            
            Shank_Angle_Stride_1_Hilbert = hilbert(Shank_Angle_Stride_1_Centered);
            Shank_Angle_Stride_2_Hilbert = hilbert(Shank_Angle_Stride_2_Centered);
            Shank_Angle_Stride_3_Hilbert = hilbert(Shank_Angle_Stride_3_Centered);
            Shank_Angle_Stride_4_Hilbert = hilbert(Shank_Angle_Stride_4_Centered);
            Shank_Angle_Stride_5_Hilbert = hilbert(Shank_Angle_Stride_5_Centered);
           

                %Imaginary part extraction
            Foot_Angle_Stride_1_Hilbert_Imag = imag(Foot_Angle_Stride_1_Hilbert);
            Foot_Angle_Stride_2_Hilbert_Imag = imag(Foot_Angle_Stride_2_Hilbert);
            Foot_Angle_Stride_3_Hilbert_Imag = imag(Foot_Angle_Stride_3_Hilbert);
            Foot_Angle_Stride_4_Hilbert_Imag = imag(Foot_Angle_Stride_4_Hilbert);
            Foot_Angle_Stride_5_Hilbert_Imag = imag(Foot_Angle_Stride_5_Hilbert);

            Shank_Angle_Stride_1_Hilbert_Imag = imag(Shank_Angle_Stride_1_Hilbert);
            Shank_Angle_Stride_2_Hilbert_Imag = imag(Shank_Angle_Stride_2_Hilbert);
            Shank_Angle_Stride_3_Hilbert_Imag = imag(Shank_Angle_Stride_3_Hilbert);
            Shank_Angle_Stride_4_Hilbert_Imag = imag(Shank_Angle_Stride_4_Hilbert);
            Shank_Angle_Stride_5_Hilbert_Imag = imag(Shank_Angle_Stride_5_Hilbert);


            % Step 3 Phase angle calculation
            for i = 1: length(Foot_Angle_Stride_1)
                Phase_Angle_Foot_Stride_1(i) = atan (Foot_Angle_Stride_1_Hilbert_Imag(i) ./ Foot_Angle_Stride_1_Centered(i) );
                Phase_Angle_Foot_Stride_2(i) = atan (Foot_Angle_Stride_2_Hilbert_Imag(i) ./ Foot_Angle_Stride_2_Centered(i) );
                Phase_Angle_Foot_Stride_3(i) = atan (Foot_Angle_Stride_3_Hilbert_Imag(i) ./ Foot_Angle_Stride_3_Centered(i) );
                Phase_Angle_Foot_Stride_4(i) = atan (Foot_Angle_Stride_4_Hilbert_Imag(i) ./ Foot_Angle_Stride_4_Centered(i) );
                Phase_Angle_Foot_Stride_5(i) = atan (Foot_Angle_Stride_5_Hilbert_Imag(i) ./ Foot_Angle_Stride_5_Centered(i) );
                                
                Phase_Angle_Shank_Stride_1(i) = atan (Shank_Angle_Stride_1_Hilbert_Imag(i) ./ Shank_Angle_Stride_1_Centered(i) );
                Phase_Angle_Shank_Stride_2(i) = atan (Shank_Angle_Stride_2_Hilbert_Imag(i) ./ Shank_Angle_Stride_2_Centered(i) );
                Phase_Angle_Shank_Stride_3(i) = atan (Shank_Angle_Stride_3_Hilbert_Imag(i) ./ Shank_Angle_Stride_3_Centered(i) );
                Phase_Angle_Shank_Stride_4(i) = atan (Shank_Angle_Stride_4_Hilbert_Imag(i) ./ Shank_Angle_Stride_4_Centered(i) );
                Phase_Angle_Shank_Stride_5(i) = atan (Shank_Angle_Stride_5_Hilbert_Imag(i) ./ Shank_Angle_Stride_5_Centered(i) );
                
            end

            % Step 4 Continous Relative Phase (CRP) calculation (Proximal-Distal)
            for i = 1: length(Foot_Angle_Stride_1)
                CRP_ShankFoot_Stride_1(i) = Phase_Angle_Shank_Stride_1(i) - Phase_Angle_Foot_Stride_1(i);
                CRP_ShankFoot_Stride_2(i) = Phase_Angle_Shank_Stride_2(i) - Phase_Angle_Foot_Stride_2(i);
                CRP_ShankFoot_Stride_3(i) = Phase_Angle_Shank_Stride_3(i) - Phase_Angle_Foot_Stride_3(i);
                CRP_ShankFoot_Stride_4(i) = Phase_Angle_Shank_Stride_4(i) - Phase_Angle_Foot_Stride_4(i);
                CRP_ShankFoot_Stride_5(i) = Phase_Angle_Shank_Stride_5(i) - Phase_Angle_Foot_Stride_5(i);
            end
            
            CRP_ShankFoot_Stride_1 = CRP_ShankFoot_Stride_1';
            CRP_ShankFoot_Stride_2 = CRP_ShankFoot_Stride_2';
            CRP_ShankFoot_Stride_3 = CRP_ShankFoot_Stride_3';
            CRP_ShankFoot_Stride_4 = CRP_ShankFoot_Stride_4';
            CRP_ShankFoot_Stride_5 = CRP_ShankFoot_Stride_5';
            
            CRP_ShankFoot_Stride_1_Mean = mean(CRP_ShankFoot_Stride_1);
            CRP_ShankFoot_Stride_2_Mean = mean(CRP_ShankFoot_Stride_2);
            CRP_ShankFoot_Stride_3_Mean = mean(CRP_ShankFoot_Stride_3);
            CRP_ShankFoot_Stride_4_Mean = mean(CRP_ShankFoot_Stride_4);
            CRP_ShankFoot_Stride_5_Mean = mean(CRP_ShankFoot_Stride_5);
            
            CRP_ShankFoot_AllStrides_Matrix_Mean = [CRP_ShankFoot_Stride_1_Mean, CRP_ShankFoot_Stride_2_Mean, CRP_ShankFoot_Stride_3_Mean, CRP_ShankFoot_Stride_4_Mean, CRP_ShankFoot_Stride_5_Mean];
            CRP_ShankFoot_Total_Mean = mean(CRP_ShankFoot_AllStrides_Matrix_Mean);
                        
            CRP_ShankFoot_Stride_1_SD = std(CRP_ShankFoot_Stride_1);
            CRP_ShankFoot_Stride_2_SD = std(CRP_ShankFoot_Stride_2);
            CRP_ShankFoot_Stride_3_SD = std(CRP_ShankFoot_Stride_3);
            CRP_ShankFoot_Stride_4_SD = std(CRP_ShankFoot_Stride_4);
            CRP_ShankFoot_Stride_5_SD = std(CRP_ShankFoot_Stride_5);
            
            CRP_ShankFoot_AllStrides_Matrix_SD = [CRP_ShankFoot_Stride_1_SD, CRP_ShankFoot_Stride_2_SD, CRP_ShankFoot_Stride_3_SD, CRP_ShankFoot_Stride_4_SD, CRP_ShankFoot_Stride_5_SD];
            CRP_ShankFoot_Total_SD = mean(CRP_ShankFoot_AllStrides_Matrix_SD);
    
            
               %One row from the FINAL export matrix per input file
               
                  %CRP Curves from each stride
                  CRP_Curves_ShankFoot_Export =[CRP_ShankFoot_Stride_1, CRP_ShankFoot_Stride_2, CRP_ShankFoot_Stride_3, CRP_ShankFoot_Stride_4, CRP_ShankFoot_Stride_5];
   
                    %Empty matrix.  I am creating this empty matrix so I can
                    % creat one single file to Export in the end for all
                    % depths conditions. The depths that does not have the
                    % segments/joints data will have empty (NaN) cells
                    
                  CRP_Curves_ThighShank_Export = [NaN(100,1), NaN(100,1), NaN(100,1), NaN(100,1), NaN(100,1)];
                  CRP_Curves_TrunkThigh_Export = [NaN(100,1), NaN(100,1), NaN(100,1), NaN(100,1), NaN(100,1)];
                  
                  CRP_Curves_KneeAnkle_Export = [NaN(100,1), NaN(100,1), NaN(100,1), NaN(100,1), NaN(100,1)];
                  CRP_Curves_HipKnee_Export = [NaN(100,1), NaN(100,1), NaN(100,1), NaN(100,1), NaN(100,1)];
                  
                  
                  %Calculating the Mean along the 5 strides from the CRP
                    %curves 0-100% to export
                  CRP_Curves_ShankFoot_Export_Mean(:,k) = mean(CRP_Curves_ShankFoot_Export,2);
                  CRP_Curves_ThighShank_Export_Mean(:,k) = mean(CRP_Curves_ThighShank_Export,2);
                  CRP_Curves_TrunkThigh_Export_Mean(:,k) = mean(CRP_Curves_TrunkThigh_Export,2);
                  
                  CRP_Curves_KneeAnkle_Export_Mean(:,k) = mean(CRP_Curves_KneeAnkle_Export,2);
                  CRP_Curves_HipKnee_Export_Mean(:,k) = mean(CRP_Curves_HipKnee_Export,2);
                   
                   %CRP Mean and SD value
                  CRP_ShankFoot_Export = [CRP_ShankFoot_Total_Mean, CRP_ShankFoot_Total_SD];
                  
                    %Empty matrix. To export in the end one single file for
                    %all depths.
                  CRP_ThighShank_Export = [nan, nan];
                  CRP_TrunkThigh_Export = [nan, nan];

                  CRP_KneeAnkle_Export = [nan, nan];
                  CRP_HipKnee_Export = [nan, nan];
                                    
                  
                  CRP_Matrix_Export(k,:) = [k, CRP_ShankFoot_Export, CRP_ThighShank_Export, CRP_TrunkThigh_Export, CRP_KneeAnkle_Export, CRP_HipKnee_Export];
       
                   
    elseif Depth == 2

        
            % Hip Depth
        Foot_Angle = Data(:,1:5);
        Shank_Angle = Data(:,7:11);
        Thigh_Angle = Data(:,13:17);

        Ankle_Joint_Angle = Data(:,19:23);
        Knee_Joint_Angle = Data(:,25:29);

                %Separating per Strides
        Foot_Angle_Stride_1 = Foot_Angle(:,1);
        Foot_Angle_Stride_2 = Foot_Angle(:,2);
        Foot_Angle_Stride_3 = Foot_Angle(:,3);
        Foot_Angle_Stride_4 = Foot_Angle(:,4);
        Foot_Angle_Stride_5 = Foot_Angle(:,5);

        Shank_Angle_Stride_1 = Shank_Angle(:,1);
        Shank_Angle_Stride_2 = Shank_Angle(:,2);
        Shank_Angle_Stride_3 = Shank_Angle(:,3);
        Shank_Angle_Stride_4 = Shank_Angle(:,4);
        Shank_Angle_Stride_5 = Shank_Angle(:,5);

        Thigh_Angle_Stride_1 = Thigh_Angle(:,1);
        Thigh_Angle_Stride_2 = Thigh_Angle(:,2);
        Thigh_Angle_Stride_3 = Thigh_Angle(:,3);
        Thigh_Angle_Stride_4 = Thigh_Angle(:,4);
        Thigh_Angle_Stride_5 = Thigh_Angle(:,5);

        Ankle_Joint_Angle_Stride_1 = Ankle_Joint_Angle(:,1);
        Ankle_Joint_Angle_Stride_2 = Ankle_Joint_Angle(:,2);
        Ankle_Joint_Angle_Stride_3 = Ankle_Joint_Angle(:,3);
        Ankle_Joint_Angle_Stride_4 = Ankle_Joint_Angle(:,4);
        Ankle_Joint_Angle_Stride_5 = Ankle_Joint_Angle(:,5);

        Knee_Joint_Angle_Stride_1 = Knee_Joint_Angle(:,1);
        Knee_Joint_Angle_Stride_2 = Knee_Joint_Angle(:,2);
        Knee_Joint_Angle_Stride_3 = Knee_Joint_Angle(:,3);
        Knee_Joint_Angle_Stride_4 = Knee_Joint_Angle(:,4);
        Knee_Joint_Angle_Stride_5 = Knee_Joint_Angle(:,5);

        
                % CRP Calculation

            % Step 1 Centering angle
            for i = 1:length(Foot_Angle_Stride_1)
                Foot_Angle_Stride_1_Centered(i) =  Foot_Angle_Stride_1(i) - min(Foot_Angle_Stride_1) -  (( max(Foot_Angle_Stride_1) - min(Foot_Angle_Stride_1) )  /2);
                Foot_Angle_Stride_2_Centered(i) =  Foot_Angle_Stride_2(i) - min(Foot_Angle_Stride_2) -  (( max(Foot_Angle_Stride_2) - min(Foot_Angle_Stride_2) )  /2);
                Foot_Angle_Stride_3_Centered(i) =  Foot_Angle_Stride_3(i) - min(Foot_Angle_Stride_3) -  (( max(Foot_Angle_Stride_3) - min(Foot_Angle_Stride_3) )  /2);
                Foot_Angle_Stride_4_Centered(i) =  Foot_Angle_Stride_4(i) - min(Foot_Angle_Stride_4) -  (( max(Foot_Angle_Stride_4) - min(Foot_Angle_Stride_4) )  /2);
                Foot_Angle_Stride_5_Centered(i) =  Foot_Angle_Stride_5(i) - min(Foot_Angle_Stride_5) -  (( max(Foot_Angle_Stride_5) - min(Foot_Angle_Stride_5) )  /2);
                
                Shank_Angle_Stride_1_Centered(i) =  Shank_Angle_Stride_1(i) - min(Shank_Angle_Stride_1) -  (( max(Shank_Angle_Stride_1) - min(Shank_Angle_Stride_1) )  /2);
                Shank_Angle_Stride_2_Centered(i) =  Shank_Angle_Stride_2(i) - min(Shank_Angle_Stride_2) -  (( max(Shank_Angle_Stride_2) - min(Shank_Angle_Stride_2) )  /2);
                Shank_Angle_Stride_3_Centered(i) =  Shank_Angle_Stride_3(i) - min(Shank_Angle_Stride_3) -  (( max(Shank_Angle_Stride_3) - min(Shank_Angle_Stride_3) )  /2);
                Shank_Angle_Stride_4_Centered(i) =  Shank_Angle_Stride_4(i) - min(Shank_Angle_Stride_4) -  (( max(Shank_Angle_Stride_4) - min(Shank_Angle_Stride_4) )  /2);
                Shank_Angle_Stride_5_Centered(i) =  Shank_Angle_Stride_5(i) - min(Shank_Angle_Stride_5) -  (( max(Shank_Angle_Stride_5) - min(Shank_Angle_Stride_5) )  /2);
                
                Thigh_Angle_Stride_1_Centered(i) =  Thigh_Angle_Stride_1(i) - min(Thigh_Angle_Stride_1) -  (( max(Thigh_Angle_Stride_1) - min(Thigh_Angle_Stride_1) )  /2);
                Thigh_Angle_Stride_2_Centered(i) =  Thigh_Angle_Stride_2(i) - min(Thigh_Angle_Stride_2) -  (( max(Thigh_Angle_Stride_2) - min(Thigh_Angle_Stride_2) )  /2);
                Thigh_Angle_Stride_3_Centered(i) =  Thigh_Angle_Stride_3(i) - min(Thigh_Angle_Stride_3) -  (( max(Thigh_Angle_Stride_3) - min(Thigh_Angle_Stride_3) )  /2);
                Thigh_Angle_Stride_4_Centered(i) =  Thigh_Angle_Stride_4(i) - min(Thigh_Angle_Stride_4) -  (( max(Thigh_Angle_Stride_4) - min(Thigh_Angle_Stride_4) )  /2);
                Thigh_Angle_Stride_5_Centered(i) =  Thigh_Angle_Stride_5(i) - min(Thigh_Angle_Stride_5) -  (( max(Thigh_Angle_Stride_5) - min(Thigh_Angle_Stride_5) )  /2);

                Ankle_Joint_Angle_Stride_1_Centered(i) =  Ankle_Joint_Angle_Stride_1(i) - min(Ankle_Joint_Angle_Stride_1) -  (( max(Ankle_Joint_Angle_Stride_1) - min(Ankle_Joint_Angle_Stride_1) )  /2);
                Ankle_Joint_Angle_Stride_2_Centered(i) =  Ankle_Joint_Angle_Stride_2(i) - min(Ankle_Joint_Angle_Stride_2) -  (( max(Ankle_Joint_Angle_Stride_2) - min(Ankle_Joint_Angle_Stride_2) )  /2);
                Ankle_Joint_Angle_Stride_3_Centered(i) =  Ankle_Joint_Angle_Stride_3(i) - min(Ankle_Joint_Angle_Stride_3) -  (( max(Ankle_Joint_Angle_Stride_3) - min(Ankle_Joint_Angle_Stride_3) )  /2);
                Ankle_Joint_Angle_Stride_4_Centered(i) =  Ankle_Joint_Angle_Stride_4(i) - min(Ankle_Joint_Angle_Stride_4) -  (( max(Ankle_Joint_Angle_Stride_4) - min(Ankle_Joint_Angle_Stride_4) )  /2);
                Ankle_Joint_Angle_Stride_5_Centered(i) =  Ankle_Joint_Angle_Stride_5(i) - min(Ankle_Joint_Angle_Stride_5) -  (( max(Ankle_Joint_Angle_Stride_5) - min(Ankle_Joint_Angle_Stride_5) )  /2);
                
                Knee_Joint_Angle_Stride_1_Centered(i) =  Knee_Joint_Angle_Stride_1(i) - min(Knee_Joint_Angle_Stride_1) -  (( max(Knee_Joint_Angle_Stride_1) - min(Knee_Joint_Angle_Stride_1) )  /2);
                Knee_Joint_Angle_Stride_2_Centered(i) =  Knee_Joint_Angle_Stride_2(i) - min(Knee_Joint_Angle_Stride_2) -  (( max(Knee_Joint_Angle_Stride_2) - min(Knee_Joint_Angle_Stride_2) )  /2);
                Knee_Joint_Angle_Stride_3_Centered(i) =  Knee_Joint_Angle_Stride_3(i) - min(Knee_Joint_Angle_Stride_3) -  (( max(Knee_Joint_Angle_Stride_3) - min(Knee_Joint_Angle_Stride_3) )  /2);
                Knee_Joint_Angle_Stride_4_Centered(i) =  Knee_Joint_Angle_Stride_4(i) - min(Knee_Joint_Angle_Stride_4) -  (( max(Knee_Joint_Angle_Stride_4) - min(Knee_Joint_Angle_Stride_4) )  /2);
                Knee_Joint_Angle_Stride_5_Centered(i) =  Knee_Joint_Angle_Stride_5(i) - min(Knee_Joint_Angle_Stride_5) -  (( max(Knee_Joint_Angle_Stride_5) - min(Knee_Joint_Angle_Stride_5) )  /2);

                
            end

            % Step 2 Hilbert transform
            Foot_Angle_Stride_1_Hilbert = hilbert(Foot_Angle_Stride_1_Centered);
            Foot_Angle_Stride_2_Hilbert = hilbert(Foot_Angle_Stride_2_Centered);
            Foot_Angle_Stride_3_Hilbert = hilbert(Foot_Angle_Stride_3_Centered);
            Foot_Angle_Stride_4_Hilbert = hilbert(Foot_Angle_Stride_4_Centered);
            Foot_Angle_Stride_5_Hilbert = hilbert(Foot_Angle_Stride_5_Centered);
            
            Shank_Angle_Stride_1_Hilbert = hilbert(Shank_Angle_Stride_1_Centered);
            Shank_Angle_Stride_2_Hilbert = hilbert(Shank_Angle_Stride_2_Centered);
            Shank_Angle_Stride_3_Hilbert = hilbert(Shank_Angle_Stride_3_Centered);
            Shank_Angle_Stride_4_Hilbert = hilbert(Shank_Angle_Stride_4_Centered);
            Shank_Angle_Stride_5_Hilbert = hilbert(Shank_Angle_Stride_5_Centered);
            
            Thigh_Angle_Stride_1_Hilbert = hilbert(Thigh_Angle_Stride_1_Centered);
            Thigh_Angle_Stride_2_Hilbert = hilbert(Thigh_Angle_Stride_2_Centered);
            Thigh_Angle_Stride_3_Hilbert = hilbert(Thigh_Angle_Stride_3_Centered);
            Thigh_Angle_Stride_4_Hilbert = hilbert(Thigh_Angle_Stride_4_Centered);
            Thigh_Angle_Stride_5_Hilbert = hilbert(Thigh_Angle_Stride_5_Centered);
           
            Ankle_Joint_Angle_Stride_1_Hilbert = hilbert(Ankle_Joint_Angle_Stride_1_Centered);
            Ankle_Joint_Angle_Stride_2_Hilbert = hilbert(Ankle_Joint_Angle_Stride_2_Centered);
            Ankle_Joint_Angle_Stride_3_Hilbert = hilbert(Ankle_Joint_Angle_Stride_3_Centered);
            Ankle_Joint_Angle_Stride_4_Hilbert = hilbert(Ankle_Joint_Angle_Stride_4_Centered);
            Ankle_Joint_Angle_Stride_5_Hilbert = hilbert(Ankle_Joint_Angle_Stride_5_Centered);
                        
            Knee_Joint_Angle_Stride_1_Hilbert = hilbert(Knee_Joint_Angle_Stride_1_Centered);
            Knee_Joint_Angle_Stride_2_Hilbert = hilbert(Knee_Joint_Angle_Stride_2_Centered);
            Knee_Joint_Angle_Stride_3_Hilbert = hilbert(Knee_Joint_Angle_Stride_3_Centered);
            Knee_Joint_Angle_Stride_4_Hilbert = hilbert(Knee_Joint_Angle_Stride_4_Centered);
            Knee_Joint_Angle_Stride_5_Hilbert = hilbert(Knee_Joint_Angle_Stride_5_Centered);
           

                %Imaginary part extraction
            Foot_Angle_Stride_1_Hilbert_Imag = imag(Foot_Angle_Stride_1_Hilbert);
            Foot_Angle_Stride_2_Hilbert_Imag = imag(Foot_Angle_Stride_2_Hilbert);
            Foot_Angle_Stride_3_Hilbert_Imag = imag(Foot_Angle_Stride_3_Hilbert);
            Foot_Angle_Stride_4_Hilbert_Imag = imag(Foot_Angle_Stride_4_Hilbert);
            Foot_Angle_Stride_5_Hilbert_Imag = imag(Foot_Angle_Stride_5_Hilbert);

            Shank_Angle_Stride_1_Hilbert_Imag = imag(Shank_Angle_Stride_1_Hilbert);
            Shank_Angle_Stride_2_Hilbert_Imag = imag(Shank_Angle_Stride_2_Hilbert);
            Shank_Angle_Stride_3_Hilbert_Imag = imag(Shank_Angle_Stride_3_Hilbert);
            Shank_Angle_Stride_4_Hilbert_Imag = imag(Shank_Angle_Stride_4_Hilbert);
            Shank_Angle_Stride_5_Hilbert_Imag = imag(Shank_Angle_Stride_5_Hilbert);
            
            Thigh_Angle_Stride_1_Hilbert_Imag = imag(Thigh_Angle_Stride_1_Hilbert);
            Thigh_Angle_Stride_2_Hilbert_Imag = imag(Thigh_Angle_Stride_2_Hilbert);
            Thigh_Angle_Stride_3_Hilbert_Imag = imag(Thigh_Angle_Stride_3_Hilbert);
            Thigh_Angle_Stride_4_Hilbert_Imag = imag(Thigh_Angle_Stride_4_Hilbert);
            Thigh_Angle_Stride_5_Hilbert_Imag = imag(Thigh_Angle_Stride_5_Hilbert);
            
            Ankle_Joint_Angle_Stride_1_Hilbert_Imag = imag(Ankle_Joint_Angle_Stride_1_Hilbert);
            Ankle_Joint_Angle_Stride_2_Hilbert_Imag = imag(Ankle_Joint_Angle_Stride_2_Hilbert);
            Ankle_Joint_Angle_Stride_3_Hilbert_Imag = imag(Ankle_Joint_Angle_Stride_3_Hilbert);
            Ankle_Joint_Angle_Stride_4_Hilbert_Imag = imag(Ankle_Joint_Angle_Stride_4_Hilbert);
            Ankle_Joint_Angle_Stride_5_Hilbert_Imag = imag(Ankle_Joint_Angle_Stride_5_Hilbert);
            
            Knee_Joint_Angle_Stride_1_Hilbert_Imag = imag(Knee_Joint_Angle_Stride_1_Hilbert);
            Knee_Joint_Angle_Stride_2_Hilbert_Imag = imag(Knee_Joint_Angle_Stride_2_Hilbert);
            Knee_Joint_Angle_Stride_3_Hilbert_Imag = imag(Knee_Joint_Angle_Stride_3_Hilbert);
            Knee_Joint_Angle_Stride_4_Hilbert_Imag = imag(Knee_Joint_Angle_Stride_4_Hilbert);
            Knee_Joint_Angle_Stride_5_Hilbert_Imag = imag(Knee_Joint_Angle_Stride_5_Hilbert);


            % Step 3 Phase angle calculation
            for i = 1: length(Foot_Angle_Stride_1)
                Phase_Angle_Foot_Stride_1(i) = atan (Foot_Angle_Stride_1_Hilbert_Imag(i) ./ Foot_Angle_Stride_1_Centered(i) );
                Phase_Angle_Foot_Stride_2(i) = atan (Foot_Angle_Stride_2_Hilbert_Imag(i) ./ Foot_Angle_Stride_2_Centered(i) );
                Phase_Angle_Foot_Stride_3(i) = atan (Foot_Angle_Stride_3_Hilbert_Imag(i) ./ Foot_Angle_Stride_3_Centered(i) );
                Phase_Angle_Foot_Stride_4(i) = atan (Foot_Angle_Stride_4_Hilbert_Imag(i) ./ Foot_Angle_Stride_4_Centered(i) );
                Phase_Angle_Foot_Stride_5(i) = atan (Foot_Angle_Stride_5_Hilbert_Imag(i) ./ Foot_Angle_Stride_5_Centered(i) );
                                
                Phase_Angle_Shank_Stride_1(i) = atan (Shank_Angle_Stride_1_Hilbert_Imag(i) ./ Shank_Angle_Stride_1_Centered(i) );
                Phase_Angle_Shank_Stride_2(i) = atan (Shank_Angle_Stride_2_Hilbert_Imag(i) ./ Shank_Angle_Stride_2_Centered(i) );
                Phase_Angle_Shank_Stride_3(i) = atan (Shank_Angle_Stride_3_Hilbert_Imag(i) ./ Shank_Angle_Stride_3_Centered(i) );
                Phase_Angle_Shank_Stride_4(i) = atan (Shank_Angle_Stride_4_Hilbert_Imag(i) ./ Shank_Angle_Stride_4_Centered(i) );
                Phase_Angle_Shank_Stride_5(i) = atan (Shank_Angle_Stride_5_Hilbert_Imag(i) ./ Shank_Angle_Stride_5_Centered(i) );
                
                Phase_Angle_Thigh_Stride_1(i) = atan (Thigh_Angle_Stride_1_Hilbert_Imag(i) ./ Thigh_Angle_Stride_1_Centered(i) );
                Phase_Angle_Thigh_Stride_2(i) = atan (Thigh_Angle_Stride_2_Hilbert_Imag(i) ./ Thigh_Angle_Stride_2_Centered(i) );
                Phase_Angle_Thigh_Stride_3(i) = atan (Thigh_Angle_Stride_3_Hilbert_Imag(i) ./ Thigh_Angle_Stride_3_Centered(i) );
                Phase_Angle_Thigh_Stride_4(i) = atan (Thigh_Angle_Stride_4_Hilbert_Imag(i) ./ Thigh_Angle_Stride_4_Centered(i) );
                Phase_Angle_Thigh_Stride_5(i) = atan (Thigh_Angle_Stride_5_Hilbert_Imag(i) ./ Thigh_Angle_Stride_5_Centered(i) );
                
                Phase_Angle_Ankle_Joint_Stride_1(i) = atan (Ankle_Joint_Angle_Stride_1_Hilbert_Imag(i) ./ Ankle_Joint_Angle_Stride_1_Centered(i) );
                Phase_Angle_Ankle_Joint_Stride_2(i) = atan (Ankle_Joint_Angle_Stride_2_Hilbert_Imag(i) ./ Ankle_Joint_Angle_Stride_2_Centered(i) );
                Phase_Angle_Ankle_Joint_Stride_3(i) = atan (Ankle_Joint_Angle_Stride_3_Hilbert_Imag(i) ./ Ankle_Joint_Angle_Stride_3_Centered(i) );
                Phase_Angle_Ankle_Joint_Stride_4(i) = atan (Ankle_Joint_Angle_Stride_4_Hilbert_Imag(i) ./ Ankle_Joint_Angle_Stride_4_Centered(i) );
                Phase_Angle_Ankle_Joint_Stride_5(i) = atan (Ankle_Joint_Angle_Stride_5_Hilbert_Imag(i) ./ Ankle_Joint_Angle_Stride_5_Centered(i) );
                
                Phase_Angle_Knee_Joint_Stride_1(i) = atan (Knee_Joint_Angle_Stride_1_Hilbert_Imag(i) ./ Knee_Joint_Angle_Stride_1_Centered(i) );
                Phase_Angle_Knee_Joint_Stride_2(i) = atan (Knee_Joint_Angle_Stride_2_Hilbert_Imag(i) ./ Knee_Joint_Angle_Stride_2_Centered(i) );
                Phase_Angle_Knee_Joint_Stride_3(i) = atan (Knee_Joint_Angle_Stride_3_Hilbert_Imag(i) ./ Knee_Joint_Angle_Stride_3_Centered(i) );
                Phase_Angle_Knee_Joint_Stride_4(i) = atan (Knee_Joint_Angle_Stride_4_Hilbert_Imag(i) ./ Knee_Joint_Angle_Stride_4_Centered(i) );
                Phase_Angle_Knee_Joint_Stride_5(i) = atan (Knee_Joint_Angle_Stride_5_Hilbert_Imag(i) ./ Knee_Joint_Angle_Stride_5_Centered(i) );
                
            end

            % Step 4 Continous Relative Phase (CRP) calculation (Proximal-Distal)
            for i = 1: length(Foot_Angle_Stride_1)
                CRP_ShankFoot_Stride_1(i) = Phase_Angle_Shank_Stride_1(i) - Phase_Angle_Foot_Stride_1(i);
                CRP_ShankFoot_Stride_2(i) = Phase_Angle_Shank_Stride_2(i) - Phase_Angle_Foot_Stride_2(i);
                CRP_ShankFoot_Stride_3(i) = Phase_Angle_Shank_Stride_3(i) - Phase_Angle_Foot_Stride_3(i);
                CRP_ShankFoot_Stride_4(i) = Phase_Angle_Shank_Stride_4(i) - Phase_Angle_Foot_Stride_4(i);
                CRP_ShankFoot_Stride_5(i) = Phase_Angle_Shank_Stride_5(i) - Phase_Angle_Foot_Stride_5(i);
                
                CRP_ThighShank_Stride_1(i) = Phase_Angle_Thigh_Stride_1(i) - Phase_Angle_Shank_Stride_1(i);
                CRP_ThighShank_Stride_2(i) = Phase_Angle_Thigh_Stride_2(i) - Phase_Angle_Shank_Stride_2(i);
                CRP_ThighShank_Stride_3(i) = Phase_Angle_Thigh_Stride_3(i) - Phase_Angle_Shank_Stride_3(i);
                CRP_ThighShank_Stride_4(i) = Phase_Angle_Thigh_Stride_4(i) - Phase_Angle_Shank_Stride_4(i);
                CRP_ThighShank_Stride_5(i) = Phase_Angle_Thigh_Stride_5(i) - Phase_Angle_Shank_Stride_5(i);
                
                CRP_KneeAnkle_Stride_1(i) = Phase_Angle_Knee_Joint_Stride_1(i) - Phase_Angle_Ankle_Joint_Stride_1(i);
                CRP_KneeAnkle_Stride_2(i) = Phase_Angle_Knee_Joint_Stride_2(i) - Phase_Angle_Ankle_Joint_Stride_2(i);
                CRP_KneeAnkle_Stride_3(i) = Phase_Angle_Knee_Joint_Stride_3(i) - Phase_Angle_Ankle_Joint_Stride_3(i);
                CRP_KneeAnkle_Stride_4(i) = Phase_Angle_Knee_Joint_Stride_4(i) - Phase_Angle_Ankle_Joint_Stride_4(i);
                CRP_KneeAnkle_Stride_5(i) = Phase_Angle_Knee_Joint_Stride_5(i) - Phase_Angle_Ankle_Joint_Stride_5(i);                
                
            end

            CRP_ShankFoot_Stride_1 = CRP_ShankFoot_Stride_1';
            CRP_ShankFoot_Stride_2 = CRP_ShankFoot_Stride_2';
            CRP_ShankFoot_Stride_3 = CRP_ShankFoot_Stride_3';
            CRP_ShankFoot_Stride_4 = CRP_ShankFoot_Stride_4';
            CRP_ShankFoot_Stride_5 = CRP_ShankFoot_Stride_5';
                        
            CRP_ThighShank_Stride_1 = CRP_ThighShank_Stride_1';
            CRP_ThighShank_Stride_2 = CRP_ThighShank_Stride_2';
            CRP_ThighShank_Stride_3 = CRP_ThighShank_Stride_3';
            CRP_ThighShank_Stride_4 = CRP_ThighShank_Stride_4';
            CRP_ThighShank_Stride_5 = CRP_ThighShank_Stride_5';
            
            CRP_KneeAnkle_Stride_1 = CRP_KneeAnkle_Stride_1';
            CRP_KneeAnkle_Stride_2 = CRP_KneeAnkle_Stride_2';
            CRP_KneeAnkle_Stride_3 = CRP_KneeAnkle_Stride_3';
            CRP_KneeAnkle_Stride_4 = CRP_KneeAnkle_Stride_4';
            CRP_KneeAnkle_Stride_5 = CRP_KneeAnkle_Stride_5';
            
            
            % Mean and SD values for each CRP pair
            
            %ShankFoot
            CRP_ShankFoot_Stride_1_Mean = mean(CRP_ShankFoot_Stride_1);
            CRP_ShankFoot_Stride_2_Mean = mean(CRP_ShankFoot_Stride_2);
            CRP_ShankFoot_Stride_3_Mean = mean(CRP_ShankFoot_Stride_3);
            CRP_ShankFoot_Stride_4_Mean = mean(CRP_ShankFoot_Stride_4);
            CRP_ShankFoot_Stride_5_Mean = mean(CRP_ShankFoot_Stride_5);
            
            CRP_ShankFoot_AllStrides_Matrix_Mean = [CRP_ShankFoot_Stride_1_Mean, CRP_ShankFoot_Stride_2_Mean, CRP_ShankFoot_Stride_3_Mean, CRP_ShankFoot_Stride_4_Mean, CRP_ShankFoot_Stride_5_Mean];
            CRP_ShankFoot_Total_Mean = mean(CRP_ShankFoot_AllStrides_Matrix_Mean);
                        
            CRP_ShankFoot_Stride_1_SD = std(CRP_ShankFoot_Stride_1);
            CRP_ShankFoot_Stride_2_SD = std(CRP_ShankFoot_Stride_2);
            CRP_ShankFoot_Stride_3_SD = std(CRP_ShankFoot_Stride_3);
            CRP_ShankFoot_Stride_4_SD = std(CRP_ShankFoot_Stride_4);
            CRP_ShankFoot_Stride_5_SD = std(CRP_ShankFoot_Stride_5);
            
            CRP_ShankFoot_AllStrides_Matrix_SD = [CRP_ShankFoot_Stride_1_SD, CRP_ShankFoot_Stride_2_SD, CRP_ShankFoot_Stride_3_SD, CRP_ShankFoot_Stride_4_SD, CRP_ShankFoot_Stride_5_SD];
            CRP_ShankFoot_Total_SD = mean(CRP_ShankFoot_AllStrides_Matrix_SD);
            
            
            %ThighShank
            CRP_ThighShank_Stride_1_Mean = mean(CRP_ThighShank_Stride_1);
            CRP_ThighShank_Stride_2_Mean = mean(CRP_ThighShank_Stride_2);
            CRP_ThighShank_Stride_3_Mean = mean(CRP_ThighShank_Stride_3);
            CRP_ThighShank_Stride_4_Mean = mean(CRP_ThighShank_Stride_4);
            CRP_ThighShank_Stride_5_Mean = mean(CRP_ThighShank_Stride_5);
            
            CRP_ThighShank_AllStrides_Matrix_Mean = [CRP_ThighShank_Stride_1_Mean, CRP_ThighShank_Stride_2_Mean, CRP_ThighShank_Stride_3_Mean, CRP_ThighShank_Stride_4_Mean, CRP_ThighShank_Stride_5_Mean];
            CRP_ThighShank_Total_Mean = mean(CRP_ThighShank_AllStrides_Matrix_Mean);
                        
            CRP_ThighShank_Stride_1_SD = std(CRP_ThighShank_Stride_1);
            CRP_ThighShank_Stride_2_SD = std(CRP_ThighShank_Stride_2);
            CRP_ThighShank_Stride_3_SD = std(CRP_ThighShank_Stride_3);
            CRP_ThighShank_Stride_4_SD = std(CRP_ThighShank_Stride_4);
            CRP_ThighShank_Stride_5_SD = std(CRP_ThighShank_Stride_5);
            
            CRP_ThighShank_AllStrides_Matrix_SD = [CRP_ThighShank_Stride_1_SD, CRP_ThighShank_Stride_2_SD, CRP_ThighShank_Stride_3_SD, CRP_ThighShank_Stride_4_SD, CRP_ThighShank_Stride_5_SD];
            CRP_ThighShank_Total_SD = mean(CRP_ThighShank_AllStrides_Matrix_SD);
    
            %KneeAnkle
            CRP_KneeAnkle_Stride_1_Mean = mean(CRP_KneeAnkle_Stride_1);
            CRP_KneeAnkle_Stride_2_Mean = mean(CRP_KneeAnkle_Stride_2);
            CRP_KneeAnkle_Stride_3_Mean = mean(CRP_KneeAnkle_Stride_3);
            CRP_KneeAnkle_Stride_4_Mean = mean(CRP_KneeAnkle_Stride_4);
            CRP_KneeAnkle_Stride_5_Mean = mean(CRP_KneeAnkle_Stride_5);
            
            CRP_KneeAnkle_AllStrides_Matrix_Mean = [CRP_KneeAnkle_Stride_1_Mean, CRP_KneeAnkle_Stride_2_Mean, CRP_KneeAnkle_Stride_3_Mean, CRP_KneeAnkle_Stride_4_Mean, CRP_KneeAnkle_Stride_5_Mean];
            CRP_KneeAnkle_Total_Mean = mean(CRP_KneeAnkle_AllStrides_Matrix_Mean);
                        
            CRP_KneeAnkle_Stride_1_SD = std(CRP_KneeAnkle_Stride_1);
            CRP_KneeAnkle_Stride_2_SD = std(CRP_KneeAnkle_Stride_2);
            CRP_KneeAnkle_Stride_3_SD = std(CRP_KneeAnkle_Stride_3);
            CRP_KneeAnkle_Stride_4_SD = std(CRP_KneeAnkle_Stride_4);
            CRP_KneeAnkle_Stride_5_SD = std(CRP_KneeAnkle_Stride_5);
            
            CRP_KneeAnkle_AllStrides_Matrix_SD = [CRP_KneeAnkle_Stride_1_SD, CRP_KneeAnkle_Stride_2_SD, CRP_KneeAnkle_Stride_3_SD, CRP_KneeAnkle_Stride_4_SD, CRP_KneeAnkle_Stride_5_SD];
            CRP_KneeAnkle_Total_SD = mean(CRP_KneeAnkle_AllStrides_Matrix_SD);
    
            
            
               %Prepare matrix to export: (EDIT HERE THE DEPTHS THAT DO NOT
               %HAVE ALL SEGMENTS OR JOINTS)
               
                  %CRP Curves from each stride
                  CRP_Curves_ShankFoot_Export =[CRP_ShankFoot_Stride_1, CRP_ShankFoot_Stride_2, CRP_ShankFoot_Stride_3, CRP_ShankFoot_Stride_4, CRP_ShankFoot_Stride_5];
                  CRP_Curves_ThighShank_Export =[CRP_ThighShank_Stride_1, CRP_ThighShank_Stride_2, CRP_ThighShank_Stride_3, CRP_ThighShank_Stride_4, CRP_ThighShank_Stride_5];
                  CRP_Curves_KneeAnkle_Export =[CRP_KneeAnkle_Stride_1, CRP_KneeAnkle_Stride_2, CRP_KneeAnkle_Stride_3, CRP_KneeAnkle_Stride_4, CRP_KneeAnkle_Stride_5];
                  
                  %Blank matrix 
                  CRP_Curves_TrunkThigh_Export = [NaN(100,1), NaN(100,1), NaN(100,1), NaN(100,1), NaN(100,1)];
                  CRP_Curves_HipKnee_Export = [NaN(100,1), NaN(100,1), NaN(100,1), NaN(100,1), NaN(100,1)];
                
                    %Calculating the Mean along the 5 strides from the CRP
                    %curves 0-100% to export
                  CRP_Curves_ShankFoot_Export_Mean(:,k) = mean(CRP_Curves_ShankFoot_Export,2);
                  CRP_Curves_ThighShank_Export_Mean(:,k) = mean(CRP_Curves_ThighShank_Export,2);
                  CRP_Curves_TrunkThigh_Export_Mean(:,k) = mean(CRP_Curves_TrunkThigh_Export,2);
                  
                  CRP_Curves_KneeAnkle_Export_Mean(:,k) = mean(CRP_Curves_KneeAnkle_Export,2);
                  CRP_Curves_HipKnee_Export_Mean(:,k) = mean(CRP_Curves_HipKnee_Export,2);
                                    
                   
                   %CRP Mean and SD value
                   CRP_ShankFoot_Export = [CRP_ShankFoot_Total_Mean, CRP_ShankFoot_Total_SD];
                   CRP_ThighShank_Export = [CRP_ThighShank_Total_Mean, CRP_ThighShank_Total_SD];
                   CRP_KneeAnkle_Export = [CRP_KneeAnkle_Total_Mean, CRP_KneeAnkle_Total_SD];
                   
                   %Blank matrix
                   CRP_TrunkThigh_Export = [nan, nan];
                   CRP_HipKnee_Export = [nan, nan];
  
                   CRP_Matrix_Export(k,:) = [k, CRP_ShankFoot_Export, CRP_ThighShank_Export, CRP_TrunkThigh_Export, CRP_KneeAnkle_Export, CRP_HipKnee_Export];
                   
    
    elseif Depth == 3 || Depth == 4

        % Umbilical & Xiphoid Depths
        Foot_Angle = Data(:,1:5);
        Shank_Angle = Data(:,7:11);
        Thigh_Angle = Data(:,13:17);
        Trunk_Angle = Data(:,19:23);

        Ankle_Joint_Angle = Data(:,25:29);
        Knee_Joint_Angle = Data(:,31:35);
        Hip_Joint_Angle = Data(:,37:41);

                %Separating per Strides
        Foot_Angle_Stride_1 = Foot_Angle(:,1);
        Foot_Angle_Stride_2 = Foot_Angle(:,2);
        Foot_Angle_Stride_3 = Foot_Angle(:,3);
        Foot_Angle_Stride_4 = Foot_Angle(:,4);
        Foot_Angle_Stride_5 = Foot_Angle(:,5);

        Shank_Angle_Stride_1 = Shank_Angle(:,1);
        Shank_Angle_Stride_2 = Shank_Angle(:,2);
        Shank_Angle_Stride_3 = Shank_Angle(:,3);
        Shank_Angle_Stride_4 = Shank_Angle(:,4);
        Shank_Angle_Stride_5 = Shank_Angle(:,5);

        Thigh_Angle_Stride_1 = Thigh_Angle(:,1);
        Thigh_Angle_Stride_2 = Thigh_Angle(:,2);
        Thigh_Angle_Stride_3 = Thigh_Angle(:,3);
        Thigh_Angle_Stride_4 = Thigh_Angle(:,4);
        Thigh_Angle_Stride_5 = Thigh_Angle(:,5);

        Trunk_Angle_Stride_1 = Trunk_Angle(:,1);
        Trunk_Angle_Stride_2 = Trunk_Angle(:,2);
        Trunk_Angle_Stride_3 = Trunk_Angle(:,3);
        Trunk_Angle_Stride_4 = Trunk_Angle(:,4);
        Trunk_Angle_Stride_5 = Trunk_Angle(:,5);


        Ankle_Joint_Angle_Stride_1 = Ankle_Joint_Angle(:,1);
        Ankle_Joint_Angle_Stride_2 = Ankle_Joint_Angle(:,2);
        Ankle_Joint_Angle_Stride_3 = Ankle_Joint_Angle(:,3);
        Ankle_Joint_Angle_Stride_4 = Ankle_Joint_Angle(:,4);
        Ankle_Joint_Angle_Stride_5 = Ankle_Joint_Angle(:,5);

        Knee_Joint_Angle_Stride_1 = Knee_Joint_Angle(:,1);
        Knee_Joint_Angle_Stride_2 = Knee_Joint_Angle(:,2);
        Knee_Joint_Angle_Stride_3 = Knee_Joint_Angle(:,3);
        Knee_Joint_Angle_Stride_4 = Knee_Joint_Angle(:,4);
        Knee_Joint_Angle_Stride_5 = Knee_Joint_Angle(:,5);

        Hip_Joint_Angle_Stride_1 = Hip_Joint_Angle(:,1);
        Hip_Joint_Angle_Stride_2 = Hip_Joint_Angle(:,2);
        Hip_Joint_Angle_Stride_3 = Hip_Joint_Angle(:,3);
        Hip_Joint_Angle_Stride_4 = Hip_Joint_Angle(:,4);
        Hip_Joint_Angle_Stride_5 = Hip_Joint_Angle(:,5);

   


                    % CRP Calculation

            % Step 1 Centering angle
            for i = 1:length(Foot_Angle_Stride_1)
                Foot_Angle_Stride_1_Centered(i) =  Foot_Angle_Stride_1(i) - min(Foot_Angle_Stride_1) -  (( max(Foot_Angle_Stride_1) - min(Foot_Angle_Stride_1) )  /2);
                Foot_Angle_Stride_2_Centered(i) =  Foot_Angle_Stride_2(i) - min(Foot_Angle_Stride_2) -  (( max(Foot_Angle_Stride_2) - min(Foot_Angle_Stride_2) )  /2);
                Foot_Angle_Stride_3_Centered(i) =  Foot_Angle_Stride_3(i) - min(Foot_Angle_Stride_3) -  (( max(Foot_Angle_Stride_3) - min(Foot_Angle_Stride_3) )  /2);
                Foot_Angle_Stride_4_Centered(i) =  Foot_Angle_Stride_4(i) - min(Foot_Angle_Stride_4) -  (( max(Foot_Angle_Stride_4) - min(Foot_Angle_Stride_4) )  /2);
                Foot_Angle_Stride_5_Centered(i) =  Foot_Angle_Stride_5(i) - min(Foot_Angle_Stride_5) -  (( max(Foot_Angle_Stride_5) - min(Foot_Angle_Stride_5) )  /2);
                
                Shank_Angle_Stride_1_Centered(i) =  Shank_Angle_Stride_1(i) - min(Shank_Angle_Stride_1) -  (( max(Shank_Angle_Stride_1) - min(Shank_Angle_Stride_1) )  /2);
                Shank_Angle_Stride_2_Centered(i) =  Shank_Angle_Stride_2(i) - min(Shank_Angle_Stride_2) -  (( max(Shank_Angle_Stride_2) - min(Shank_Angle_Stride_2) )  /2);
                Shank_Angle_Stride_3_Centered(i) =  Shank_Angle_Stride_3(i) - min(Shank_Angle_Stride_3) -  (( max(Shank_Angle_Stride_3) - min(Shank_Angle_Stride_3) )  /2);
                Shank_Angle_Stride_4_Centered(i) =  Shank_Angle_Stride_4(i) - min(Shank_Angle_Stride_4) -  (( max(Shank_Angle_Stride_4) - min(Shank_Angle_Stride_4) )  /2);
                Shank_Angle_Stride_5_Centered(i) =  Shank_Angle_Stride_5(i) - min(Shank_Angle_Stride_5) -  (( max(Shank_Angle_Stride_5) - min(Shank_Angle_Stride_5) )  /2);
                
                Thigh_Angle_Stride_1_Centered(i) =  Thigh_Angle_Stride_1(i) - min(Thigh_Angle_Stride_1) -  (( max(Thigh_Angle_Stride_1) - min(Thigh_Angle_Stride_1) )  /2);
                Thigh_Angle_Stride_2_Centered(i) =  Thigh_Angle_Stride_2(i) - min(Thigh_Angle_Stride_2) -  (( max(Thigh_Angle_Stride_2) - min(Thigh_Angle_Stride_2) )  /2);
                Thigh_Angle_Stride_3_Centered(i) =  Thigh_Angle_Stride_3(i) - min(Thigh_Angle_Stride_3) -  (( max(Thigh_Angle_Stride_3) - min(Thigh_Angle_Stride_3) )  /2);
                Thigh_Angle_Stride_4_Centered(i) =  Thigh_Angle_Stride_4(i) - min(Thigh_Angle_Stride_4) -  (( max(Thigh_Angle_Stride_4) - min(Thigh_Angle_Stride_4) )  /2);
                Thigh_Angle_Stride_5_Centered(i) =  Thigh_Angle_Stride_5(i) - min(Thigh_Angle_Stride_5) -  (( max(Thigh_Angle_Stride_5) - min(Thigh_Angle_Stride_5) )  /2);

                Trunk_Angle_Stride_1_Centered(i) =  Trunk_Angle_Stride_1(i) - min(Trunk_Angle_Stride_1) -  (( max(Trunk_Angle_Stride_1) - min(Trunk_Angle_Stride_1) )  /2);
                Trunk_Angle_Stride_2_Centered(i) =  Trunk_Angle_Stride_2(i) - min(Trunk_Angle_Stride_2) -  (( max(Trunk_Angle_Stride_2) - min(Trunk_Angle_Stride_2) )  /2);
                Trunk_Angle_Stride_3_Centered(i) =  Trunk_Angle_Stride_3(i) - min(Trunk_Angle_Stride_3) -  (( max(Trunk_Angle_Stride_3) - min(Trunk_Angle_Stride_3) )  /2);
                Trunk_Angle_Stride_4_Centered(i) =  Trunk_Angle_Stride_4(i) - min(Trunk_Angle_Stride_4) -  (( max(Trunk_Angle_Stride_4) - min(Trunk_Angle_Stride_4) )  /2);
                Trunk_Angle_Stride_5_Centered(i) =  Trunk_Angle_Stride_5(i) - min(Trunk_Angle_Stride_5) -  (( max(Trunk_Angle_Stride_5) - min(Trunk_Angle_Stride_5) )  /2);

                
                Ankle_Joint_Angle_Stride_1_Centered(i) =  Ankle_Joint_Angle_Stride_1(i) - min(Ankle_Joint_Angle_Stride_1) -  (( max(Ankle_Joint_Angle_Stride_1) - min(Ankle_Joint_Angle_Stride_1) )  /2);
                Ankle_Joint_Angle_Stride_2_Centered(i) =  Ankle_Joint_Angle_Stride_2(i) - min(Ankle_Joint_Angle_Stride_2) -  (( max(Ankle_Joint_Angle_Stride_2) - min(Ankle_Joint_Angle_Stride_2) )  /2);
                Ankle_Joint_Angle_Stride_3_Centered(i) =  Ankle_Joint_Angle_Stride_3(i) - min(Ankle_Joint_Angle_Stride_3) -  (( max(Ankle_Joint_Angle_Stride_3) - min(Ankle_Joint_Angle_Stride_3) )  /2);
                Ankle_Joint_Angle_Stride_4_Centered(i) =  Ankle_Joint_Angle_Stride_4(i) - min(Ankle_Joint_Angle_Stride_4) -  (( max(Ankle_Joint_Angle_Stride_4) - min(Ankle_Joint_Angle_Stride_4) )  /2);
                Ankle_Joint_Angle_Stride_5_Centered(i) =  Ankle_Joint_Angle_Stride_5(i) - min(Ankle_Joint_Angle_Stride_5) -  (( max(Ankle_Joint_Angle_Stride_5) - min(Ankle_Joint_Angle_Stride_5) )  /2);
                
                Knee_Joint_Angle_Stride_1_Centered(i) =  Knee_Joint_Angle_Stride_1(i) - min(Knee_Joint_Angle_Stride_1) -  (( max(Knee_Joint_Angle_Stride_1) - min(Knee_Joint_Angle_Stride_1) )  /2);
                Knee_Joint_Angle_Stride_2_Centered(i) =  Knee_Joint_Angle_Stride_2(i) - min(Knee_Joint_Angle_Stride_2) -  (( max(Knee_Joint_Angle_Stride_2) - min(Knee_Joint_Angle_Stride_2) )  /2);
                Knee_Joint_Angle_Stride_3_Centered(i) =  Knee_Joint_Angle_Stride_3(i) - min(Knee_Joint_Angle_Stride_3) -  (( max(Knee_Joint_Angle_Stride_3) - min(Knee_Joint_Angle_Stride_3) )  /2);
                Knee_Joint_Angle_Stride_4_Centered(i) =  Knee_Joint_Angle_Stride_4(i) - min(Knee_Joint_Angle_Stride_4) -  (( max(Knee_Joint_Angle_Stride_4) - min(Knee_Joint_Angle_Stride_4) )  /2);
                Knee_Joint_Angle_Stride_5_Centered(i) =  Knee_Joint_Angle_Stride_5(i) - min(Knee_Joint_Angle_Stride_5) -  (( max(Knee_Joint_Angle_Stride_5) - min(Knee_Joint_Angle_Stride_5) )  /2);
                
                Hip_Joint_Angle_Stride_1_Centered(i) =  Hip_Joint_Angle_Stride_1(i) - min(Hip_Joint_Angle_Stride_1) -  (( max(Hip_Joint_Angle_Stride_1) - min(Hip_Joint_Angle_Stride_1) )  /2);
                Hip_Joint_Angle_Stride_2_Centered(i) =  Hip_Joint_Angle_Stride_2(i) - min(Hip_Joint_Angle_Stride_2) -  (( max(Hip_Joint_Angle_Stride_2) - min(Hip_Joint_Angle_Stride_2) )  /2);
                Hip_Joint_Angle_Stride_3_Centered(i) =  Hip_Joint_Angle_Stride_3(i) - min(Hip_Joint_Angle_Stride_3) -  (( max(Hip_Joint_Angle_Stride_3) - min(Hip_Joint_Angle_Stride_3) )  /2);
                Hip_Joint_Angle_Stride_4_Centered(i) =  Hip_Joint_Angle_Stride_4(i) - min(Hip_Joint_Angle_Stride_4) -  (( max(Hip_Joint_Angle_Stride_4) - min(Hip_Joint_Angle_Stride_4) )  /2);
                Hip_Joint_Angle_Stride_5_Centered(i) =  Hip_Joint_Angle_Stride_5(i) - min(Hip_Joint_Angle_Stride_5) -  (( max(Hip_Joint_Angle_Stride_5) - min(Hip_Joint_Angle_Stride_5) )  /2);

            end

            % Step 2 Hilbert transform
            Foot_Angle_Stride_1_Hilbert = hilbert(Foot_Angle_Stride_1_Centered);
            Foot_Angle_Stride_2_Hilbert = hilbert(Foot_Angle_Stride_2_Centered);
            Foot_Angle_Stride_3_Hilbert = hilbert(Foot_Angle_Stride_3_Centered);
            Foot_Angle_Stride_4_Hilbert = hilbert(Foot_Angle_Stride_4_Centered);
            Foot_Angle_Stride_5_Hilbert = hilbert(Foot_Angle_Stride_5_Centered);
            
            Shank_Angle_Stride_1_Hilbert = hilbert(Shank_Angle_Stride_1_Centered);
            Shank_Angle_Stride_2_Hilbert = hilbert(Shank_Angle_Stride_2_Centered);
            Shank_Angle_Stride_3_Hilbert = hilbert(Shank_Angle_Stride_3_Centered);
            Shank_Angle_Stride_4_Hilbert = hilbert(Shank_Angle_Stride_4_Centered);
            Shank_Angle_Stride_5_Hilbert = hilbert(Shank_Angle_Stride_5_Centered);
            
            Thigh_Angle_Stride_1_Hilbert = hilbert(Thigh_Angle_Stride_1_Centered);
            Thigh_Angle_Stride_2_Hilbert = hilbert(Thigh_Angle_Stride_2_Centered);
            Thigh_Angle_Stride_3_Hilbert = hilbert(Thigh_Angle_Stride_3_Centered);
            Thigh_Angle_Stride_4_Hilbert = hilbert(Thigh_Angle_Stride_4_Centered);
            Thigh_Angle_Stride_5_Hilbert = hilbert(Thigh_Angle_Stride_5_Centered);
            
            Trunk_Angle_Stride_1_Hilbert = hilbert(Trunk_Angle_Stride_1_Centered);
            Trunk_Angle_Stride_2_Hilbert = hilbert(Trunk_Angle_Stride_2_Centered);
            Trunk_Angle_Stride_3_Hilbert = hilbert(Trunk_Angle_Stride_3_Centered);
            Trunk_Angle_Stride_4_Hilbert = hilbert(Trunk_Angle_Stride_4_Centered);
            Trunk_Angle_Stride_5_Hilbert = hilbert(Trunk_Angle_Stride_5_Centered);
           
            Ankle_Joint_Angle_Stride_1_Hilbert = hilbert(Ankle_Joint_Angle_Stride_1_Centered);
            Ankle_Joint_Angle_Stride_2_Hilbert = hilbert(Ankle_Joint_Angle_Stride_2_Centered);
            Ankle_Joint_Angle_Stride_3_Hilbert = hilbert(Ankle_Joint_Angle_Stride_3_Centered);
            Ankle_Joint_Angle_Stride_4_Hilbert = hilbert(Ankle_Joint_Angle_Stride_4_Centered);
            Ankle_Joint_Angle_Stride_5_Hilbert = hilbert(Ankle_Joint_Angle_Stride_5_Centered);
                        
            Knee_Joint_Angle_Stride_1_Hilbert = hilbert(Knee_Joint_Angle_Stride_1_Centered);
            Knee_Joint_Angle_Stride_2_Hilbert = hilbert(Knee_Joint_Angle_Stride_2_Centered);
            Knee_Joint_Angle_Stride_3_Hilbert = hilbert(Knee_Joint_Angle_Stride_3_Centered);
            Knee_Joint_Angle_Stride_4_Hilbert = hilbert(Knee_Joint_Angle_Stride_4_Centered);
            Knee_Joint_Angle_Stride_5_Hilbert = hilbert(Knee_Joint_Angle_Stride_5_Centered);
            
            Hip_Joint_Angle_Stride_1_Hilbert = hilbert(Hip_Joint_Angle_Stride_1_Centered);
            Hip_Joint_Angle_Stride_2_Hilbert = hilbert(Hip_Joint_Angle_Stride_2_Centered);
            Hip_Joint_Angle_Stride_3_Hilbert = hilbert(Hip_Joint_Angle_Stride_3_Centered);
            Hip_Joint_Angle_Stride_4_Hilbert = hilbert(Hip_Joint_Angle_Stride_4_Centered);
            Hip_Joint_Angle_Stride_5_Hilbert = hilbert(Hip_Joint_Angle_Stride_5_Centered);
           

                %Imaginary part extraction
            Foot_Angle_Stride_1_Hilbert_Imag = imag(Foot_Angle_Stride_1_Hilbert);
            Foot_Angle_Stride_2_Hilbert_Imag = imag(Foot_Angle_Stride_2_Hilbert);
            Foot_Angle_Stride_3_Hilbert_Imag = imag(Foot_Angle_Stride_3_Hilbert);
            Foot_Angle_Stride_4_Hilbert_Imag = imag(Foot_Angle_Stride_4_Hilbert);
            Foot_Angle_Stride_5_Hilbert_Imag = imag(Foot_Angle_Stride_5_Hilbert);

            Shank_Angle_Stride_1_Hilbert_Imag = imag(Shank_Angle_Stride_1_Hilbert);
            Shank_Angle_Stride_2_Hilbert_Imag = imag(Shank_Angle_Stride_2_Hilbert);
            Shank_Angle_Stride_3_Hilbert_Imag = imag(Shank_Angle_Stride_3_Hilbert);
            Shank_Angle_Stride_4_Hilbert_Imag = imag(Shank_Angle_Stride_4_Hilbert);
            Shank_Angle_Stride_5_Hilbert_Imag = imag(Shank_Angle_Stride_5_Hilbert);
            
            Thigh_Angle_Stride_1_Hilbert_Imag = imag(Thigh_Angle_Stride_1_Hilbert);
            Thigh_Angle_Stride_2_Hilbert_Imag = imag(Thigh_Angle_Stride_2_Hilbert);
            Thigh_Angle_Stride_3_Hilbert_Imag = imag(Thigh_Angle_Stride_3_Hilbert);
            Thigh_Angle_Stride_4_Hilbert_Imag = imag(Thigh_Angle_Stride_4_Hilbert);
            Thigh_Angle_Stride_5_Hilbert_Imag = imag(Thigh_Angle_Stride_5_Hilbert);
            
            Trunk_Angle_Stride_1_Hilbert_Imag = imag(Trunk_Angle_Stride_1_Hilbert);
            Trunk_Angle_Stride_2_Hilbert_Imag = imag(Trunk_Angle_Stride_2_Hilbert);
            Trunk_Angle_Stride_3_Hilbert_Imag = imag(Trunk_Angle_Stride_3_Hilbert);
            Trunk_Angle_Stride_4_Hilbert_Imag = imag(Trunk_Angle_Stride_4_Hilbert);
            Trunk_Angle_Stride_5_Hilbert_Imag = imag(Trunk_Angle_Stride_5_Hilbert);
            
            Ankle_Joint_Angle_Stride_1_Hilbert_Imag = imag(Ankle_Joint_Angle_Stride_1_Hilbert);
            Ankle_Joint_Angle_Stride_2_Hilbert_Imag = imag(Ankle_Joint_Angle_Stride_2_Hilbert);
            Ankle_Joint_Angle_Stride_3_Hilbert_Imag = imag(Ankle_Joint_Angle_Stride_3_Hilbert);
            Ankle_Joint_Angle_Stride_4_Hilbert_Imag = imag(Ankle_Joint_Angle_Stride_4_Hilbert);
            Ankle_Joint_Angle_Stride_5_Hilbert_Imag = imag(Ankle_Joint_Angle_Stride_5_Hilbert);
            
            Knee_Joint_Angle_Stride_1_Hilbert_Imag = imag(Knee_Joint_Angle_Stride_1_Hilbert);
            Knee_Joint_Angle_Stride_2_Hilbert_Imag = imag(Knee_Joint_Angle_Stride_2_Hilbert);
            Knee_Joint_Angle_Stride_3_Hilbert_Imag = imag(Knee_Joint_Angle_Stride_3_Hilbert);
            Knee_Joint_Angle_Stride_4_Hilbert_Imag = imag(Knee_Joint_Angle_Stride_4_Hilbert);
            Knee_Joint_Angle_Stride_5_Hilbert_Imag = imag(Knee_Joint_Angle_Stride_5_Hilbert);
            
            Hip_Joint_Angle_Stride_1_Hilbert_Imag = imag(Hip_Joint_Angle_Stride_1_Hilbert);
            Hip_Joint_Angle_Stride_2_Hilbert_Imag = imag(Hip_Joint_Angle_Stride_2_Hilbert);
            Hip_Joint_Angle_Stride_3_Hilbert_Imag = imag(Hip_Joint_Angle_Stride_3_Hilbert);
            Hip_Joint_Angle_Stride_4_Hilbert_Imag = imag(Hip_Joint_Angle_Stride_4_Hilbert);
            Hip_Joint_Angle_Stride_5_Hilbert_Imag = imag(Hip_Joint_Angle_Stride_5_Hilbert);


            % Step 3 Phase angle calculation
            for i = 1: length(Foot_Angle_Stride_1)
                Phase_Angle_Foot_Stride_1(i) = atan (Foot_Angle_Stride_1_Hilbert_Imag(i) ./ Foot_Angle_Stride_1_Centered(i) );
                Phase_Angle_Foot_Stride_2(i) = atan (Foot_Angle_Stride_2_Hilbert_Imag(i) ./ Foot_Angle_Stride_2_Centered(i) );
                Phase_Angle_Foot_Stride_3(i) = atan (Foot_Angle_Stride_3_Hilbert_Imag(i) ./ Foot_Angle_Stride_3_Centered(i) );
                Phase_Angle_Foot_Stride_4(i) = atan (Foot_Angle_Stride_4_Hilbert_Imag(i) ./ Foot_Angle_Stride_4_Centered(i) );
                Phase_Angle_Foot_Stride_5(i) = atan (Foot_Angle_Stride_5_Hilbert_Imag(i) ./ Foot_Angle_Stride_5_Centered(i) );
                                
                Phase_Angle_Shank_Stride_1(i) = atan (Shank_Angle_Stride_1_Hilbert_Imag(i) ./ Shank_Angle_Stride_1_Centered(i) );
                Phase_Angle_Shank_Stride_2(i) = atan (Shank_Angle_Stride_2_Hilbert_Imag(i) ./ Shank_Angle_Stride_2_Centered(i) );
                Phase_Angle_Shank_Stride_3(i) = atan (Shank_Angle_Stride_3_Hilbert_Imag(i) ./ Shank_Angle_Stride_3_Centered(i) );
                Phase_Angle_Shank_Stride_4(i) = atan (Shank_Angle_Stride_4_Hilbert_Imag(i) ./ Shank_Angle_Stride_4_Centered(i) );
                Phase_Angle_Shank_Stride_5(i) = atan (Shank_Angle_Stride_5_Hilbert_Imag(i) ./ Shank_Angle_Stride_5_Centered(i) );
                
                Phase_Angle_Thigh_Stride_1(i) = atan (Thigh_Angle_Stride_1_Hilbert_Imag(i) ./ Thigh_Angle_Stride_1_Centered(i) );
                Phase_Angle_Thigh_Stride_2(i) = atan (Thigh_Angle_Stride_2_Hilbert_Imag(i) ./ Thigh_Angle_Stride_2_Centered(i) );
                Phase_Angle_Thigh_Stride_3(i) = atan (Thigh_Angle_Stride_3_Hilbert_Imag(i) ./ Thigh_Angle_Stride_3_Centered(i) );
                Phase_Angle_Thigh_Stride_4(i) = atan (Thigh_Angle_Stride_4_Hilbert_Imag(i) ./ Thigh_Angle_Stride_4_Centered(i) );
                Phase_Angle_Thigh_Stride_5(i) = atan (Thigh_Angle_Stride_5_Hilbert_Imag(i) ./ Thigh_Angle_Stride_5_Centered(i) );
                
                Phase_Angle_Trunk_Stride_1(i) = atan (Trunk_Angle_Stride_1_Hilbert_Imag(i) ./ Trunk_Angle_Stride_1_Centered(i) );
                Phase_Angle_Trunk_Stride_2(i) = atan (Trunk_Angle_Stride_2_Hilbert_Imag(i) ./ Trunk_Angle_Stride_2_Centered(i) );
                Phase_Angle_Trunk_Stride_3(i) = atan (Trunk_Angle_Stride_3_Hilbert_Imag(i) ./ Trunk_Angle_Stride_3_Centered(i) );
                Phase_Angle_Trunk_Stride_4(i) = atan (Trunk_Angle_Stride_4_Hilbert_Imag(i) ./ Trunk_Angle_Stride_4_Centered(i) );
                Phase_Angle_Trunk_Stride_5(i) = atan (Trunk_Angle_Stride_5_Hilbert_Imag(i) ./ Trunk_Angle_Stride_5_Centered(i) );
                
                Phase_Angle_Ankle_Joint_Stride_1(i) = atan (Ankle_Joint_Angle_Stride_1_Hilbert_Imag(i) ./ Ankle_Joint_Angle_Stride_1_Centered(i) );
                Phase_Angle_Ankle_Joint_Stride_2(i) = atan (Ankle_Joint_Angle_Stride_2_Hilbert_Imag(i) ./ Ankle_Joint_Angle_Stride_2_Centered(i) );
                Phase_Angle_Ankle_Joint_Stride_3(i) = atan (Ankle_Joint_Angle_Stride_3_Hilbert_Imag(i) ./ Ankle_Joint_Angle_Stride_3_Centered(i) );
                Phase_Angle_Ankle_Joint_Stride_4(i) = atan (Ankle_Joint_Angle_Stride_4_Hilbert_Imag(i) ./ Ankle_Joint_Angle_Stride_4_Centered(i) );
                Phase_Angle_Ankle_Joint_Stride_5(i) = atan (Ankle_Joint_Angle_Stride_5_Hilbert_Imag(i) ./ Ankle_Joint_Angle_Stride_5_Centered(i) );
                
                Phase_Angle_Knee_Joint_Stride_1(i) = atan (Knee_Joint_Angle_Stride_1_Hilbert_Imag(i) ./ Knee_Joint_Angle_Stride_1_Centered(i) );
                Phase_Angle_Knee_Joint_Stride_2(i) = atan (Knee_Joint_Angle_Stride_2_Hilbert_Imag(i) ./ Knee_Joint_Angle_Stride_2_Centered(i) );
                Phase_Angle_Knee_Joint_Stride_3(i) = atan (Knee_Joint_Angle_Stride_3_Hilbert_Imag(i) ./ Knee_Joint_Angle_Stride_3_Centered(i) );
                Phase_Angle_Knee_Joint_Stride_4(i) = atan (Knee_Joint_Angle_Stride_4_Hilbert_Imag(i) ./ Knee_Joint_Angle_Stride_4_Centered(i) );
                Phase_Angle_Knee_Joint_Stride_5(i) = atan (Knee_Joint_Angle_Stride_5_Hilbert_Imag(i) ./ Knee_Joint_Angle_Stride_5_Centered(i) );
                
                Phase_Angle_Hip_Joint_Stride_1(i) = atan (Hip_Joint_Angle_Stride_1_Hilbert_Imag(i) ./ Hip_Joint_Angle_Stride_1_Centered(i) );
                Phase_Angle_Hip_Joint_Stride_2(i) = atan (Hip_Joint_Angle_Stride_2_Hilbert_Imag(i) ./ Hip_Joint_Angle_Stride_2_Centered(i) );
                Phase_Angle_Hip_Joint_Stride_3(i) = atan (Hip_Joint_Angle_Stride_3_Hilbert_Imag(i) ./ Hip_Joint_Angle_Stride_3_Centered(i) );
                Phase_Angle_Hip_Joint_Stride_4(i) = atan (Hip_Joint_Angle_Stride_4_Hilbert_Imag(i) ./ Hip_Joint_Angle_Stride_4_Centered(i) );
                Phase_Angle_Hip_Joint_Stride_5(i) = atan (Hip_Joint_Angle_Stride_5_Hilbert_Imag(i) ./ Hip_Joint_Angle_Stride_5_Centered(i) );
                
            end

            % Step 4 Continous Relative Phase (CRP) calculation (Proximal-Distal)
            for i = 1: length(Foot_Angle_Stride_1)
                CRP_ShankFoot_Stride_1(i) = Phase_Angle_Shank_Stride_1(i) - Phase_Angle_Foot_Stride_1(i);
                CRP_ShankFoot_Stride_2(i) = Phase_Angle_Shank_Stride_2(i) - Phase_Angle_Foot_Stride_2(i);
                CRP_ShankFoot_Stride_3(i) = Phase_Angle_Shank_Stride_3(i) - Phase_Angle_Foot_Stride_3(i);
                CRP_ShankFoot_Stride_4(i) = Phase_Angle_Shank_Stride_4(i) - Phase_Angle_Foot_Stride_4(i);
                CRP_ShankFoot_Stride_5(i) = Phase_Angle_Shank_Stride_5(i) - Phase_Angle_Foot_Stride_5(i);
                
                CRP_ThighShank_Stride_1(i) = Phase_Angle_Thigh_Stride_1(i) - Phase_Angle_Shank_Stride_1(i);
                CRP_ThighShank_Stride_2(i) = Phase_Angle_Thigh_Stride_2(i) - Phase_Angle_Shank_Stride_2(i);
                CRP_ThighShank_Stride_3(i) = Phase_Angle_Thigh_Stride_3(i) - Phase_Angle_Shank_Stride_3(i);
                CRP_ThighShank_Stride_4(i) = Phase_Angle_Thigh_Stride_4(i) - Phase_Angle_Shank_Stride_4(i);
                CRP_ThighShank_Stride_5(i) = Phase_Angle_Thigh_Stride_5(i) - Phase_Angle_Shank_Stride_5(i);
                
                CRP_TrunkThigh_Stride_1(i) = Phase_Angle_Trunk_Stride_1(i) - Phase_Angle_Thigh_Stride_1(i);
                CRP_TrunkThigh_Stride_2(i) = Phase_Angle_Trunk_Stride_2(i) - Phase_Angle_Thigh_Stride_2(i);
                CRP_TrunkThigh_Stride_3(i) = Phase_Angle_Trunk_Stride_3(i) - Phase_Angle_Thigh_Stride_3(i);
                CRP_TrunkThigh_Stride_4(i) = Phase_Angle_Trunk_Stride_4(i) - Phase_Angle_Thigh_Stride_4(i);
                CRP_TrunkThigh_Stride_5(i) = Phase_Angle_Trunk_Stride_5(i) - Phase_Angle_Thigh_Stride_5(i);
                
                CRP_KneeAnkle_Stride_1(i) = Phase_Angle_Knee_Joint_Stride_1(i) - Phase_Angle_Ankle_Joint_Stride_1(i);
                CRP_KneeAnkle_Stride_2(i) = Phase_Angle_Knee_Joint_Stride_2(i) - Phase_Angle_Ankle_Joint_Stride_2(i);
                CRP_KneeAnkle_Stride_3(i) = Phase_Angle_Knee_Joint_Stride_3(i) - Phase_Angle_Ankle_Joint_Stride_3(i);
                CRP_KneeAnkle_Stride_4(i) = Phase_Angle_Knee_Joint_Stride_4(i) - Phase_Angle_Ankle_Joint_Stride_4(i);
                CRP_KneeAnkle_Stride_5(i) = Phase_Angle_Knee_Joint_Stride_5(i) - Phase_Angle_Ankle_Joint_Stride_5(i);  
                
                CRP_HipKnee_Stride_1(i) = Phase_Angle_Hip_Joint_Stride_1(i) - Phase_Angle_Knee_Joint_Stride_1(i);
                CRP_HipKnee_Stride_2(i) = Phase_Angle_Hip_Joint_Stride_2(i) - Phase_Angle_Knee_Joint_Stride_2(i);
                CRP_HipKnee_Stride_3(i) = Phase_Angle_Hip_Joint_Stride_3(i) - Phase_Angle_Knee_Joint_Stride_3(i);
                CRP_HipKnee_Stride_4(i) = Phase_Angle_Hip_Joint_Stride_4(i) - Phase_Angle_Knee_Joint_Stride_4(i);
                CRP_HipKnee_Stride_5(i) = Phase_Angle_Hip_Joint_Stride_5(i) - Phase_Angle_Knee_Joint_Stride_5(i);
                
            end
            
            
            CRP_ShankFoot_Stride_1 = CRP_ShankFoot_Stride_1';
            CRP_ShankFoot_Stride_2 = CRP_ShankFoot_Stride_2';
            CRP_ShankFoot_Stride_3 = CRP_ShankFoot_Stride_3';
            CRP_ShankFoot_Stride_4 = CRP_ShankFoot_Stride_4';
            CRP_ShankFoot_Stride_5 = CRP_ShankFoot_Stride_5';
                        
            CRP_ThighShank_Stride_1 = CRP_ThighShank_Stride_1';
            CRP_ThighShank_Stride_2 = CRP_ThighShank_Stride_2';
            CRP_ThighShank_Stride_3 = CRP_ThighShank_Stride_3';
            CRP_ThighShank_Stride_4 = CRP_ThighShank_Stride_4';
            CRP_ThighShank_Stride_5 = CRP_ThighShank_Stride_5';
            
            CRP_TrunkThigh_Stride_1 = CRP_TrunkThigh_Stride_1';
            CRP_TrunkThigh_Stride_2 = CRP_TrunkThigh_Stride_2';
            CRP_TrunkThigh_Stride_3 = CRP_TrunkThigh_Stride_3';
            CRP_TrunkThigh_Stride_4 = CRP_TrunkThigh_Stride_4';
            CRP_TrunkThigh_Stride_5 = CRP_TrunkThigh_Stride_5';
            
            CRP_KneeAnkle_Stride_1 = CRP_KneeAnkle_Stride_1';
            CRP_KneeAnkle_Stride_2 = CRP_KneeAnkle_Stride_2';
            CRP_KneeAnkle_Stride_3 = CRP_KneeAnkle_Stride_3';
            CRP_KneeAnkle_Stride_4 = CRP_KneeAnkle_Stride_4';
            CRP_KneeAnkle_Stride_5 = CRP_KneeAnkle_Stride_5';
            
            CRP_HipKnee_Stride_1 = CRP_HipKnee_Stride_1';
            CRP_HipKnee_Stride_2 = CRP_HipKnee_Stride_2';
            CRP_HipKnee_Stride_3 = CRP_HipKnee_Stride_3';
            CRP_HipKnee_Stride_4 = CRP_HipKnee_Stride_4';
            CRP_HipKnee_Stride_5 = CRP_HipKnee_Stride_5';
                        

            % Mean and SD values for each CRP pair
            
            %ShankFoot
            CRP_ShankFoot_Stride_1_Mean = mean(CRP_ShankFoot_Stride_1);
            CRP_ShankFoot_Stride_2_Mean = mean(CRP_ShankFoot_Stride_2);
            CRP_ShankFoot_Stride_3_Mean = mean(CRP_ShankFoot_Stride_3);
            CRP_ShankFoot_Stride_4_Mean = mean(CRP_ShankFoot_Stride_4);
            CRP_ShankFoot_Stride_5_Mean = mean(CRP_ShankFoot_Stride_5);
            
            CRP_ShankFoot_AllStrides_Matrix_Mean = [CRP_ShankFoot_Stride_1_Mean, CRP_ShankFoot_Stride_2_Mean, CRP_ShankFoot_Stride_3_Mean, CRP_ShankFoot_Stride_4_Mean, CRP_ShankFoot_Stride_5_Mean];
            CRP_ShankFoot_Total_Mean = mean(CRP_ShankFoot_AllStrides_Matrix_Mean);
                        
            CRP_ShankFoot_Stride_1_SD = std(CRP_ShankFoot_Stride_1);
            CRP_ShankFoot_Stride_2_SD = std(CRP_ShankFoot_Stride_2);
            CRP_ShankFoot_Stride_3_SD = std(CRP_ShankFoot_Stride_3);
            CRP_ShankFoot_Stride_4_SD = std(CRP_ShankFoot_Stride_4);
            CRP_ShankFoot_Stride_5_SD = std(CRP_ShankFoot_Stride_5);
            
            CRP_ShankFoot_AllStrides_Matrix_SD = [CRP_ShankFoot_Stride_1_SD, CRP_ShankFoot_Stride_2_SD, CRP_ShankFoot_Stride_3_SD, CRP_ShankFoot_Stride_4_SD, CRP_ShankFoot_Stride_5_SD];
            CRP_ShankFoot_Total_SD = mean(CRP_ShankFoot_AllStrides_Matrix_SD);
            
            
            %ThighShank
            CRP_ThighShank_Stride_1_Mean = mean(CRP_ThighShank_Stride_1);
            CRP_ThighShank_Stride_2_Mean = mean(CRP_ThighShank_Stride_2);
            CRP_ThighShank_Stride_3_Mean = mean(CRP_ThighShank_Stride_3);
            CRP_ThighShank_Stride_4_Mean = mean(CRP_ThighShank_Stride_4);
            CRP_ThighShank_Stride_5_Mean = mean(CRP_ThighShank_Stride_5);
            
            CRP_ThighShank_AllStrides_Matrix_Mean = [CRP_ThighShank_Stride_1_Mean, CRP_ThighShank_Stride_2_Mean, CRP_ThighShank_Stride_3_Mean, CRP_ThighShank_Stride_4_Mean, CRP_ThighShank_Stride_5_Mean];
            CRP_ThighShank_Total_Mean = mean(CRP_ThighShank_AllStrides_Matrix_Mean);
                        
            CRP_ThighShank_Stride_1_SD = std(CRP_ThighShank_Stride_1);
            CRP_ThighShank_Stride_2_SD = std(CRP_ThighShank_Stride_2);
            CRP_ThighShank_Stride_3_SD = std(CRP_ThighShank_Stride_3);
            CRP_ThighShank_Stride_4_SD = std(CRP_ThighShank_Stride_4);
            CRP_ThighShank_Stride_5_SD = std(CRP_ThighShank_Stride_5);
            
            CRP_ThighShank_AllStrides_Matrix_SD = [CRP_ThighShank_Stride_1_SD, CRP_ThighShank_Stride_2_SD, CRP_ThighShank_Stride_3_SD, CRP_ThighShank_Stride_4_SD, CRP_ThighShank_Stride_5_SD];
            CRP_ThighShank_Total_SD = mean(CRP_ThighShank_AllStrides_Matrix_SD);
    
            %TrunkThigh
            CRP_TrunkThigh_Stride_1_Mean = mean(CRP_TrunkThigh_Stride_1);
            CRP_TrunkThigh_Stride_2_Mean = mean(CRP_TrunkThigh_Stride_2);
            CRP_TrunkThigh_Stride_3_Mean = mean(CRP_TrunkThigh_Stride_3);
            CRP_TrunkThigh_Stride_4_Mean = mean(CRP_TrunkThigh_Stride_4);
            CRP_TrunkThigh_Stride_5_Mean = mean(CRP_TrunkThigh_Stride_5);
            
            CRP_TrunkThigh_AllStrides_Matrix_Mean = [CRP_TrunkThigh_Stride_1_Mean, CRP_TrunkThigh_Stride_2_Mean, CRP_TrunkThigh_Stride_3_Mean, CRP_TrunkThigh_Stride_4_Mean, CRP_TrunkThigh_Stride_5_Mean];
            CRP_TrunkThigh_Total_Mean = mean(CRP_TrunkThigh_AllStrides_Matrix_Mean);
                        
            CRP_TrunkThigh_Stride_1_SD = std(CRP_TrunkThigh_Stride_1);
            CRP_TrunkThigh_Stride_2_SD = std(CRP_TrunkThigh_Stride_2);
            CRP_TrunkThigh_Stride_3_SD = std(CRP_TrunkThigh_Stride_3);
            CRP_TrunkThigh_Stride_4_SD = std(CRP_TrunkThigh_Stride_4);
            CRP_TrunkThigh_Stride_5_SD = std(CRP_TrunkThigh_Stride_5);
            
            CRP_TrunkThigh_AllStrides_Matrix_SD = [CRP_TrunkThigh_Stride_1_SD, CRP_TrunkThigh_Stride_2_SD, CRP_TrunkThigh_Stride_3_SD, CRP_TrunkThigh_Stride_4_SD, CRP_TrunkThigh_Stride_5_SD];
            CRP_TrunkThigh_Total_SD = mean(CRP_TrunkThigh_AllStrides_Matrix_SD);
    
            
            %KneeAnkle
            CRP_KneeAnkle_Stride_1_Mean = mean(CRP_KneeAnkle_Stride_1);
            CRP_KneeAnkle_Stride_2_Mean = mean(CRP_KneeAnkle_Stride_2);
            CRP_KneeAnkle_Stride_3_Mean = mean(CRP_KneeAnkle_Stride_3);
            CRP_KneeAnkle_Stride_4_Mean = mean(CRP_KneeAnkle_Stride_4);
            CRP_KneeAnkle_Stride_5_Mean = mean(CRP_KneeAnkle_Stride_5);
            
            CRP_KneeAnkle_AllStrides_Matrix_Mean = [CRP_KneeAnkle_Stride_1_Mean, CRP_KneeAnkle_Stride_2_Mean, CRP_KneeAnkle_Stride_3_Mean, CRP_KneeAnkle_Stride_4_Mean, CRP_KneeAnkle_Stride_5_Mean];
            CRP_KneeAnkle_Total_Mean = mean(CRP_KneeAnkle_AllStrides_Matrix_Mean);
                        
            CRP_KneeAnkle_Stride_1_SD = std(CRP_KneeAnkle_Stride_1);
            CRP_KneeAnkle_Stride_2_SD = std(CRP_KneeAnkle_Stride_2);
            CRP_KneeAnkle_Stride_3_SD = std(CRP_KneeAnkle_Stride_3);
            CRP_KneeAnkle_Stride_4_SD = std(CRP_KneeAnkle_Stride_4);
            CRP_KneeAnkle_Stride_5_SD = std(CRP_KneeAnkle_Stride_5);
            
            CRP_KneeAnkle_AllStrides_Matrix_SD = [CRP_KneeAnkle_Stride_1_SD, CRP_KneeAnkle_Stride_2_SD, CRP_KneeAnkle_Stride_3_SD, CRP_KneeAnkle_Stride_4_SD, CRP_KneeAnkle_Stride_5_SD];
            CRP_KneeAnkle_Total_SD = mean(CRP_KneeAnkle_AllStrides_Matrix_SD);
    
            %HipKnee
            CRP_HipKnee_Stride_1_Mean = mean(CRP_HipKnee_Stride_1);
            CRP_HipKnee_Stride_2_Mean = mean(CRP_HipKnee_Stride_2);
            CRP_HipKnee_Stride_3_Mean = mean(CRP_HipKnee_Stride_3);
            CRP_HipKnee_Stride_4_Mean = mean(CRP_HipKnee_Stride_4);
            CRP_HipKnee_Stride_5_Mean = mean(CRP_HipKnee_Stride_5);
            
            CRP_HipKnee_AllStrides_Matrix_Mean = [CRP_HipKnee_Stride_1_Mean, CRP_HipKnee_Stride_2_Mean, CRP_HipKnee_Stride_3_Mean, CRP_HipKnee_Stride_4_Mean, CRP_HipKnee_Stride_5_Mean];
            CRP_HipKnee_Total_Mean = mean(CRP_HipKnee_AllStrides_Matrix_Mean);
                        
            CRP_HipKnee_Stride_1_SD = std(CRP_HipKnee_Stride_1);
            CRP_HipKnee_Stride_2_SD = std(CRP_HipKnee_Stride_2);
            CRP_HipKnee_Stride_3_SD = std(CRP_HipKnee_Stride_3);
            CRP_HipKnee_Stride_4_SD = std(CRP_HipKnee_Stride_4);
            CRP_HipKnee_Stride_5_SD = std(CRP_HipKnee_Stride_5);
            
            CRP_HipKnee_AllStrides_Matrix_SD = [CRP_HipKnee_Stride_1_SD, CRP_HipKnee_Stride_2_SD, CRP_HipKnee_Stride_3_SD, CRP_HipKnee_Stride_4_SD, CRP_HipKnee_Stride_5_SD];
            CRP_HipKnee_Total_SD = mean(CRP_HipKnee_AllStrides_Matrix_SD);
            
            
               %Prepare matrix to export:
               
                  %CRP Curves from each stride
                  CRP_Curves_ShankFoot_Export =[CRP_ShankFoot_Stride_1, CRP_ShankFoot_Stride_2, CRP_ShankFoot_Stride_3, CRP_ShankFoot_Stride_4, CRP_ShankFoot_Stride_5];
                  CRP_Curves_ThighShank_Export =[CRP_ThighShank_Stride_1, CRP_ThighShank_Stride_2, CRP_ThighShank_Stride_3, CRP_ThighShank_Stride_4, CRP_ThighShank_Stride_5];
                  CRP_Curves_TrunkThigh_Export =[CRP_TrunkThigh_Stride_1, CRP_TrunkThigh_Stride_2, CRP_TrunkThigh_Stride_3, CRP_TrunkThigh_Stride_4, CRP_TrunkThigh_Stride_5];
                  
                  CRP_Curves_KneeAnkle_Export =[CRP_KneeAnkle_Stride_1, CRP_KneeAnkle_Stride_2, CRP_KneeAnkle_Stride_3, CRP_KneeAnkle_Stride_4, CRP_KneeAnkle_Stride_5];
                  CRP_Curves_HipKnee_Export =[CRP_HipKnee_Stride_1, CRP_HipKnee_Stride_2, CRP_HipKnee_Stride_3, CRP_HipKnee_Stride_4, CRP_HipKnee_Stride_5];
                  
                    %Calculating the Mean along the 5 strides from the CRP
                    %curves 0-100% to export
                  CRP_Curves_ShankFoot_Export_Mean(:,k) = mean(CRP_Curves_ShankFoot_Export,2);
                  CRP_Curves_ThighShank_Export_Mean(:,k) = mean(CRP_Curves_ThighShank_Export,2);
                  CRP_Curves_TrunkThigh_Export_Mean(:,k) = mean(CRP_Curves_TrunkThigh_Export,2);
                  
                  CRP_Curves_KneeAnkle_Export_Mean(:,k) = mean(CRP_Curves_KneeAnkle_Export,2);
                  CRP_Curves_HipKnee_Export_Mean(:,k) = mean(CRP_Curves_HipKnee_Export,2);
                                    
                   
                  %CRP Mean and SD value
                  CRP_ShankFoot_Export = [CRP_ShankFoot_Total_Mean, CRP_ShankFoot_Total_SD];
                  CRP_ThighShank_Export = [CRP_ThighShank_Total_Mean, CRP_ThighShank_Total_SD];
                  CRP_TrunkThigh_Export = [CRP_TrunkThigh_Total_Mean, CRP_TrunkThigh_Total_SD];

                  CRP_KneeAnkle_Export = [CRP_KneeAnkle_Total_Mean, CRP_KneeAnkle_Total_SD];
                  CRP_HipKnee_Export = [CRP_HipKnee_Total_Mean, CRP_HipKnee_Total_SD];
                                    
                  CRP_Matrix_Export(k,:) = [k, CRP_ShankFoot_Export, CRP_ThighShank_Export, CRP_TrunkThigh_Export, CRP_KneeAnkle_Export, CRP_HipKnee_Export];
       
    
      end % End from the Single File Loop

    
end % End from all Files loop


%% Export Section

% Contains:
    % Export CRP Mean and SD
    % Export CRP 0-100% Curves
    
    
    
% Export CRP Mean and SD
Output_File_Path_CRP_Mean_SD = [Output_File_Path 'CRP_Mean_SD' '.xls']; 

Export_Header_CRP_Mean_SD = {'File Number', 'Shank-Foot Mean', 'Shank Foot SD', 'Thigh-Shank Mean', 'Thigh-Shank SD', 'Trunk-Thigh Mean', 'Trunk-Thigh SD', 'Knee-Ankle Mean', 'Knee-Ankle SD', 'Hip-Knee Mean', 'Hip-Knee SD'};
CRP_Matrix_Export_Complete = [Export_Header_CRP_Mean_SD; num2cell(CRP_Matrix_Export)]; 

xlswrite(Output_File_Path_CRP_Mean_SD, CRP_Matrix_Export_Complete); 



% Export CRP 0-100% Curves

% Loop to prepare the Matrix to export CRP Curves 0-100% Files, dividing per Depth and Speed

  for k =  1:length(Input_File_Data)
    
    %Extract the depth and speed condition
    File_Name_Individual = char(Input_Matrix_FilesNames(k));
    Depth = Input_File_Data(k,2);
    Speed = Input_File_Data(k,3);
    
    if Depth == 1 && Speed == 1
        
        %Test this matrix header
        
        %Editing here%%%%%%%%%%%%%%%%%%%%%% Use this format for all others.
        %It will export as Header the number of the file. 
        CRP_Curves_ShankFoot_Export_Mean_Knee_Two(:,k) =  [(k);(CRP_Curves_ShankFoot_Export_Mean(:,k))];
        CRP_Curves_ThighShank_Export_Mean_Knee_Two(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
        CRP_Curves_TrunkThigh_Export_Mean_Knee_Two(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
        CRP_Curves_KneeAnkle_Export_Mean_Knee_Two(:,k) =  [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
        CRP_Curves_HipKnee_Export_Mean_Knee_Two(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
         
    elseif Depth == 1 && Speed == 2
         CRP_Curves_ShankFoot_Export_Mean_Knee_Four(:,k) =  [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
         CRP_Curves_ThighShank_Export_Mean_Knee_Four(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
         CRP_Curves_TrunkThigh_Export_Mean_Knee_Four(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
         CRP_Curves_KneeAnkle_Export_Mean_Knee_Four(:,k) =  [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
         CRP_Curves_HipKnee_Export_Mean_Knee_Four(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
            
    elseif Depth == 1 && Speed == 3
         CRP_Curves_ShankFoot_Export_Mean_Knee_Six(:,k) =  [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
         CRP_Curves_ThighShank_Export_Mean_Knee_Six(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
         CRP_Curves_TrunkThigh_Export_Mean_Knee_Six(:,k) =  CRP_Curves_TrunkThigh_Export_Mean(:,k);
        
         CRP_Curves_KneeAnkle_Export_Mean_Knee_Six(:,k) =  [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
         CRP_Curves_HipKnee_Export_Mean_Knee_Six(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
         
     elseif Depth == 1 && Speed == 4
         CRP_Curves_ShankFoot_Export_Mean_Knee_Eight(:,k) =  [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
         CRP_Curves_ThighShank_Export_Mean_Knee_Eight(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
         CRP_Curves_TrunkThigh_Export_Mean_Knee_Eight(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
         CRP_Curves_KneeAnkle_Export_Mean_Knee_Eight(:,k) =  [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
         CRP_Curves_HipKnee_Export_Mean_Knee_Eight(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
         
    elseif Depth == 2 && Speed == 1
        CRP_Curves_ShankFoot_Export_Mean_Hip_Two(:,k) =  [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
        CRP_Curves_ThighShank_Export_Mean_Hip_Two(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
        CRP_Curves_TrunkThigh_Export_Mean_Hip_Two(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
        CRP_Curves_KneeAnkle_Export_Mean_Hip_Two(:,k) =  [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
        CRP_Curves_HipKnee_Export_Mean_Hip_Two(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
         
    elseif Depth == 2 && Speed == 2
         CRP_Curves_ShankFoot_Export_Mean_Hip_Four(:,k) =  [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
         CRP_Curves_ThighShank_Export_Mean_Hip_Four(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
         CRP_Curves_TrunkThigh_Export_Mean_Hip_Four(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
         CRP_Curves_KneeAnkle_Export_Mean_Hip_Four(:,k) =  [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
         CRP_Curves_HipKnee_Export_Mean_Hip_Four(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
                     
    elseif Depth == 2 && Speed == 3
         CRP_Curves_ShankFoot_Export_Mean_Hip_Six(:,k) =  [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
         CRP_Curves_ThighShank_Export_Mean_Hip_Six(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
         CRP_Curves_TrunkThigh_Export_Mean_Hip_Six(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
         CRP_Curves_KneeAnkle_Export_Mean_Hip_Six(:,k) =  [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
         CRP_Curves_HipKnee_Export_Mean_Hip_Six(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
         
     elseif Depth == 2 && Speed == 4
         CRP_Curves_ShankFoot_Export_Mean_Hip_Eight(:,k) =  [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
         CRP_Curves_ThighShank_Export_Mean_Hip_Eight(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
         CRP_Curves_TrunkThigh_Export_Mean_Hip_Eight(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
         CRP_Curves_KneeAnkle_Export_Mean_Hip_Eight(:,k) =  [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
         CRP_Curves_HipKnee_Export_Mean_Hip_Eight(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
         
     elseif Depth == 3 && Speed == 1
        CRP_Curves_ShankFoot_Export_Mean_Umbilical_Two(:,k) =  [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
        CRP_Curves_ThighShank_Export_Mean_Umbilical_Two(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
        CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Two(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
        CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Two(:,k) =  [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
        CRP_Curves_HipKnee_Export_Mean_Umbilical_Two(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
         
    elseif Depth == 3 && Speed == 2
         CRP_Curves_ShankFoot_Export_Mean_Umbilical_Four(:,k) = [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
         CRP_Curves_ThighShank_Export_Mean_Umbilical_Four(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
         CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Four(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
         CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Four(:,k) = [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
         CRP_Curves_HipKnee_Export_Mean_Umbilical_Four(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
            
    elseif Depth == 3 && Speed == 3
         CRP_Curves_ShankFoot_Export_Mean_Umbilical_Six(:,k) =  [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
         CRP_Curves_ThighShank_Export_Mean_Umbilical_Six(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
         CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Six(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
         CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Six(:,k) =  [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
         CRP_Curves_HipKnee_Export_Mean_Umbilical_Six(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
         
     elseif Depth == 3 && Speed == 4
         CRP_Curves_ShankFoot_Export_Mean_Umbilical_Eight(:,k) =  [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
         CRP_Curves_ThighShank_Export_Mean_Umbilical_Eight(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
         CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Eight(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
         CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Eight(:,k) =  [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
         CRP_Curves_HipKnee_Export_Mean_Umbilical_Eight(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
              
    elseif Depth == 4 && Speed == 1
        CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Two(:,k) =  [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
        CRP_Curves_ThighShank_Export_Mean_Xiphoid_Two(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
        CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Two(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
        CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Two(:,k) = [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
        CRP_Curves_HipKnee_Export_Mean_Xiphoid_Two(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
         
    elseif Depth == 4 && Speed == 2
        CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Four(:,k) =  [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
        CRP_Curves_ThighShank_Export_Mean_Xiphoid_Four(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
        CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Four(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
        CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Four(:,k) =  [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
        CRP_Curves_HipKnee_Export_Mean_Xiphoid_Four(:,k) = [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
            
    elseif Depth == 4 && Speed == 3
         CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Six(:,k) =  [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
         CRP_Curves_ThighShank_Export_Mean_Xiphoid_Six(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
         CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Six(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
         CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Six(:,k) =  [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
         CRP_Curves_HipKnee_Export_Mean_Xiphoid_Six(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];
         
     elseif Depth == 4 && Speed == 4
         CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Eight(:,k) =  [(k);CRP_Curves_ShankFoot_Export_Mean(:,k)];
         CRP_Curves_ThighShank_Export_Mean_Xiphoid_Eight(:,k) =  [(k);CRP_Curves_ThighShank_Export_Mean(:,k)];
         CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Eight(:,k) =  [(k);CRP_Curves_TrunkThigh_Export_Mean(:,k)];
        
         CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Eight(:,k) =  [(k);CRP_Curves_KneeAnkle_Export_Mean(:,k)];
         CRP_Curves_HipKnee_Export_Mean_Xiphoid_Eight(:,k) =  [(k);CRP_Curves_HipKnee_Export_Mean(:,k)];

    end
    
    
end



%Here I delete the columns that contains only NaN and Zero values, to leave only
%the columns with actual values. But I have to test better this, to see if
%I am not deleting any column that have one single zero, but has actual
%values.
 
%Delete NaN columns
% B(:, all(isnan(B),1)) = [];

%Delete Zeros columns
% A(:,all(A == 0))=[];
   
            %     % ShankFoot
            % 
            % CRP_Curves_ShankFoot_Export_Mean_Knee_Two(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Knee_Two),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Knee_Two(:,all(CRP_Curves_ShankFoot_Export_Mean_Knee_Two == 0))=[];
            % CRP_Curves_ShankFoot_Export_Mean_Knee_Four(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Knee_Four),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Knee_Four(:,all(CRP_Curves_ShankFoot_Export_Mean_Knee_Four == 0))=[];
            % CRP_Curves_ShankFoot_Export_Mean_Knee_Six(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Knee_Six),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Knee_Six(:,all(CRP_Curves_ShankFoot_Export_Mean_Knee_Six == 0))=[];
            % CRP_Curves_ShankFoot_Export_Mean_Knee_Eight(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Knee_Eight),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Knee_Eight(:,all(CRP_Curves_ShankFoot_Export_Mean_Knee_Eight == 0))=[];
            % 
            % CRP_Curves_ShankFoot_Export_Mean_Hip_Two(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Hip_Two),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Hip_Two(:,all(CRP_Curves_ShankFoot_Export_Mean_Hip_Two == 0))=[];
            % CRP_Curves_ShankFoot_Export_Mean_Hip_Four(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Hip_Four),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Hip_Four(:,all(CRP_Curves_ShankFoot_Export_Mean_Hip_Four == 0))=[];
            % CRP_Curves_ShankFoot_Export_Mean_Hip_Six(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Hip_Six),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Hip_Six(:,all(CRP_Curves_ShankFoot_Export_Mean_Hip_Six == 0))=[];
            % CRP_Curves_ShankFoot_Export_Mean_Hip_Eight(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Hip_Eight),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Hip_Eight(:,all(CRP_Curves_ShankFoot_Export_Mean_Hip_Eight == 0))=[];
            % 
            % CRP_Curves_ShankFoot_Export_Mean_Umbilical_Two(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Umbilical_Two),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Umbilical_Two(:,all(CRP_Curves_ShankFoot_Export_Mean_Umbilical_Two == 0))=[];
            % CRP_Curves_ShankFoot_Export_Mean_Umbilical_Four(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Umbilical_Four),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Umbilical_Four(:,all(CRP_Curves_ShankFoot_Export_Mean_Umbilical_Four == 0))=[];
            % CRP_Curves_ShankFoot_Export_Mean_Umbilical_Six(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Umbilical_Six),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Umbilical_Six(:,all(CRP_Curves_ShankFoot_Export_Mean_Umbilical_Six == 0))=[];
            % CRP_Curves_ShankFoot_Export_Mean_Umbilical_Eight(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Umbilical_Eight),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Umbilical_Eight(:,all(CRP_Curves_ShankFoot_Export_Mean_Umbilical_Eight == 0))=[];
            % 
            % CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Two(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Two),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Two(:,all(CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Two == 0))=[];
            % CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Four(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Four),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Four(:,all(CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Four == 0))=[];
            % CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Six(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Six),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Six(:,all(CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Six == 0))=[];
            % CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Eight(:, all(isnan(CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Eight),1)) = [];
            % CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Eight(:,all(CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Eight == 0))=[];
            % 
            % 
            %     %ThighShank
            % CRP_Curves_ThighShank_Export_Mean_Knee_Two(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Knee_Two),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Knee_Two(:,all(CRP_Curves_ThighShank_Export_Mean_Knee_Two == 0))=[];
            % CRP_Curves_ThighShank_Export_Mean_Knee_Four(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Knee_Four),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Knee_Four(:,all(CRP_Curves_ThighShank_Export_Mean_Knee_Four == 0))=[];
            % CRP_Curves_ThighShank_Export_Mean_Knee_Six(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Knee_Six),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Knee_Six(:,all(CRP_Curves_ThighShank_Export_Mean_Knee_Six == 0))=[];
            % CRP_Curves_ThighShank_Export_Mean_Knee_Eight(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Knee_Eight),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Knee_Eight(:,all(CRP_Curves_ThighShank_Export_Mean_Knee_Eight == 0))=[];
            % 
            % CRP_Curves_ThighShank_Export_Mean_Hip_Two(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Hip_Two),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Hip_Two(:,all(CRP_Curves_ThighShank_Export_Mean_Hip_Two == 0))=[];
            % CRP_Curves_ThighShank_Export_Mean_Hip_Four(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Hip_Four),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Hip_Four(:,all(CRP_Curves_ThighShank_Export_Mean_Hip_Four == 0))=[];
            % CRP_Curves_ThighShank_Export_Mean_Hip_Six(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Hip_Six),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Hip_Six(:,all(CRP_Curves_ThighShank_Export_Mean_Hip_Six == 0))=[];
            % CRP_Curves_ThighShank_Export_Mean_Hip_Eight(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Hip_Eight),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Hip_Eight(:,all(CRP_Curves_ThighShank_Export_Mean_Hip_Eight == 0))=[];
            % 
            % CRP_Curves_ThighShank_Export_Mean_Umbilical_Two(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Umbilical_Two),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Umbilical_Two(:,all(CRP_Curves_ThighShank_Export_Mean_Umbilical_Two == 0))=[];
            % CRP_Curves_ThighShank_Export_Mean_Umbilical_Four(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Umbilical_Four),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Umbilical_Four(:,all(CRP_Curves_ThighShank_Export_Mean_Umbilical_Four == 0))=[];
            % CRP_Curves_ThighShank_Export_Mean_Umbilical_Six(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Umbilical_Six),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Umbilical_Six(:,all(CRP_Curves_ThighShank_Export_Mean_Umbilical_Six == 0))=[];
            % CRP_Curves_ThighShank_Export_Mean_Umbilical_Eight(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Umbilical_Eight),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Umbilical_Eight(:,all(CRP_Curves_ThighShank_Export_Mean_Umbilical_Eight == 0))=[];
            % 
            % CRP_Curves_ThighShank_Export_Mean_Xiphoid_Two(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Xiphoid_Two),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Xiphoid_Two(:,all(CRP_Curves_ThighShank_Export_Mean_Xiphoid_Two == 0))=[];
            % CRP_Curves_ThighShank_Export_Mean_Xiphoid_Four(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Xiphoid_Four),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Xiphoid_Four(:,all(CRP_Curves_ThighShank_Export_Mean_Xiphoid_Four == 0))=[];
            % CRP_Curves_ThighShank_Export_Mean_Xiphoid_Six(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Xiphoid_Six),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Xiphoid_Six(:,all(CRP_Curves_ThighShank_Export_Mean_Xiphoid_Six == 0))=[];
            % CRP_Curves_ThighShank_Export_Mean_Xiphoid_Eight(:, all(isnan(CRP_Curves_ThighShank_Export_Mean_Xiphoid_Eight),1)) = [];
            % CRP_Curves_ThighShank_Export_Mean_Xiphoid_Eight(:,all(CRP_Curves_ThighShank_Export_Mean_Xiphoid_Eight == 0))=[];
            % 
            %     %TrunkThigh
            % CRP_Curves_TrunkThigh_Export_Mean_Knee_Two(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Knee_Two),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Knee_Two(:,all(CRP_Curves_TrunkThigh_Export_Mean_Knee_Two == 0))=[];
            % CRP_Curves_TrunkThigh_Export_Mean_Knee_Four(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Knee_Four),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Knee_Four(:,all(CRP_Curves_TrunkThigh_Export_Mean_Knee_Four == 0))=[];
            % CRP_Curves_TrunkThigh_Export_Mean_Knee_Six(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Knee_Six),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Knee_Six(:,all(CRP_Curves_TrunkThigh_Export_Mean_Knee_Six == 0))=[];
            % CRP_Curves_TrunkThigh_Export_Mean_Knee_Eight(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Knee_Eight),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Knee_Eight(:,all(CRP_Curves_TrunkThigh_Export_Mean_Knee_Eight == 0))=[];
            % 
            % CRP_Curves_TrunkThigh_Export_Mean_Hip_Two(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Hip_Two),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Hip_Two(:,all(CRP_Curves_TrunkThigh_Export_Mean_Hip_Two == 0))=[];
            % CRP_Curves_TrunkThigh_Export_Mean_Hip_Four(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Hip_Four),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Hip_Four(:,all(CRP_Curves_TrunkThigh_Export_Mean_Hip_Four == 0))=[];
            % CRP_Curves_TrunkThigh_Export_Mean_Hip_Six(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Hip_Six),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Hip_Six(:,all(CRP_Curves_TrunkThigh_Export_Mean_Hip_Six == 0))=[];
            % CRP_Curves_TrunkThigh_Export_Mean_Hip_Eight(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Hip_Eight),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Hip_Eight(:,all(CRP_Curves_TrunkThigh_Export_Mean_Hip_Eight == 0))=[];
            % 
            % CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Two(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Two),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Two(:,all(CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Two == 0))=[];
            % CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Four(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Four),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Four(:,all(CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Four == 0))=[];
            % CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Six(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Six),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Six(:,all(CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Six == 0))=[];
            % CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Eight(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Eight),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Eight(:,all(CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Eight == 0))=[];
            % 
            % CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Two(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Two),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Two(:,all(CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Two == 0))=[];
            % CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Four(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Four),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Four(:,all(CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Four == 0))=[];
            % CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Six(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Six),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Six(:,all(CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Six == 0))=[];
            % CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Eight(:, all(isnan(CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Eight),1)) = [];
            % CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Eight(:,all(CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Eight == 0))=[];
            % 
            % 
            %     %KneeAnkle
            % CRP_Curves_KneeAnkle_Export_Mean_Knee_Two(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Knee_Two),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Knee_Two(:,all(CRP_Curves_KneeAnkle_Export_Mean_Knee_Two == 0))=[];
            % CRP_Curves_KneeAnkle_Export_Mean_Knee_Four(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Knee_Four),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Knee_Four(:,all(CRP_Curves_KneeAnkle_Export_Mean_Knee_Four == 0))=[];
            % CRP_Curves_KneeAnkle_Export_Mean_Knee_Six(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Knee_Six),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Knee_Six(:,all(CRP_Curves_KneeAnkle_Export_Mean_Knee_Six == 0))=[];
            % CRP_Curves_KneeAnkle_Export_Mean_Knee_Eight(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Knee_Eight),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Knee_Eight(:,all(CRP_Curves_KneeAnkle_Export_Mean_Knee_Eight == 0))=[];
            % 
            % CRP_Curves_KneeAnkle_Export_Mean_Hip_Two(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Hip_Two),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Hip_Two(:,all(CRP_Curves_KneeAnkle_Export_Mean_Hip_Two == 0))=[];
            % CRP_Curves_KneeAnkle_Export_Mean_Hip_Four(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Hip_Four),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Hip_Four(:,all(CRP_Curves_KneeAnkle_Export_Mean_Hip_Four == 0))=[];
            % CRP_Curves_KneeAnkle_Export_Mean_Hip_Six(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Hip_Six),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Hip_Six(:,all(CRP_Curves_KneeAnkle_Export_Mean_Hip_Six == 0))=[];
            % CRP_Curves_KneeAnkle_Export_Mean_Hip_Eight(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Hip_Eight),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Hip_Eight(:,all(CRP_Curves_KneeAnkle_Export_Mean_Hip_Eight == 0))=[];
            % 
            % CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Two(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Two),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Two(:,all(CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Two == 0))=[];
            % CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Four(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Four),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Four(:,all(CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Four == 0))=[];
            % CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Six(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Six),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Six(:,all(CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Six == 0))=[];
            % CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Eight(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Eight),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Eight(:,all(CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Eight == 0))=[];
            % 
            % CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Two(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Two),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Two(:,all(CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Two == 0))=[];
            % CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Four(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Four),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Four(:,all(CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Four == 0))=[];
            % CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Six(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Six),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Six(:,all(CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Six == 0))=[];
            % CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Eight(:, all(isnan(CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Eight),1)) = [];
            % CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Eight(:,all(CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Eight == 0))=[];
            % 
            % 
            % %HipKnee
            % 
            % CRP_Curves_HipKnee_Export_Mean_Knee_Two(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Knee_Two),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Knee_Two(:,all(CRP_Curves_HipKnee_Export_Mean_Knee_Two == 0))=[];
            % CRP_Curves_HipKnee_Export_Mean_Knee_Four(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Knee_Four),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Knee_Four(:,all(CRP_Curves_HipKnee_Export_Mean_Knee_Four == 0))=[];
            % CRP_Curves_HipKnee_Export_Mean_Knee_Six(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Knee_Six),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Knee_Six(:,all(CRP_Curves_HipKnee_Export_Mean_Knee_Six == 0))=[];
            % CRP_Curves_HipKnee_Export_Mean_Knee_Eight(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Knee_Eight),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Knee_Eight(:,all(CRP_Curves_HipKnee_Export_Mean_Knee_Eight == 0))=[];
            % 
            % CRP_Curves_HipKnee_Export_Mean_Hip_Two(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Hip_Two),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Hip_Two(:,all(CRP_Curves_HipKnee_Export_Mean_Hip_Two == 0))=[];
            % CRP_Curves_HipKnee_Export_Mean_Hip_Four(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Hip_Four),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Hip_Four(:,all(CRP_Curves_HipKnee_Export_Mean_Hip_Four == 0))=[];
            % CRP_Curves_HipKnee_Export_Mean_Hip_Six(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Hip_Six),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Hip_Six(:,all(CRP_Curves_HipKnee_Export_Mean_Hip_Six == 0))=[];
            % CRP_Curves_HipKnee_Export_Mean_Hip_Eight(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Hip_Eight),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Hip_Eight(:,all(CRP_Curves_HipKnee_Export_Mean_Hip_Eight == 0))=[];
            % 
            % CRP_Curves_HipKnee_Export_Mean_Umbilical_Two(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Umbilical_Two),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Umbilical_Two(:,all(CRP_Curves_HipKnee_Export_Mean_Umbilical_Two == 0))=[];
            % CRP_Curves_HipKnee_Export_Mean_Umbilical_Four(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Umbilical_Four),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Umbilical_Four(:,all(CRP_Curves_HipKnee_Export_Mean_Umbilical_Four == 0))=[];
            % CRP_Curves_HipKnee_Export_Mean_Umbilical_Six(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Umbilical_Six),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Umbilical_Six(:,all(CRP_Curves_HipKnee_Export_Mean_Umbilical_Six == 0))=[];
            % CRP_Curves_HipKnee_Export_Mean_Umbilical_Eight(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Umbilical_Eight),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Umbilical_Eight(:,all(CRP_Curves_HipKnee_Export_Mean_Umbilical_Eight == 0))=[];
            % 
            % CRP_Curves_HipKnee_Export_Mean_Xiphoid_Two(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Xiphoid_Two),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Xiphoid_Two(:,all(CRP_Curves_HipKnee_Export_Mean_Xiphoid_Two == 0))=[];
            % CRP_Curves_HipKnee_Export_Mean_Xiphoid_Four(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Xiphoid_Four),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Xiphoid_Four(:,all(CRP_Curves_HipKnee_Export_Mean_Xiphoid_Four == 0))=[];
            % CRP_Curves_HipKnee_Export_Mean_Xiphoid_Six(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Xiphoid_Six),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Xiphoid_Six(:,all(CRP_Curves_HipKnee_Export_Mean_Xiphoid_Six == 0))=[];
            % CRP_Curves_HipKnee_Export_Mean_Xiphoid_Eight(:, all(isnan(CRP_Curves_HipKnee_Export_Mean_Xiphoid_Eight),1)) = [];
            % CRP_Curves_HipKnee_Export_Mean_Xiphoid_Eight(:,all(CRP_Curves_HipKnee_Export_Mean_Xiphoid_Eight == 0))=[];


%Export CRP Curves

    %Export CRP Curves ShankFoot
Output_File_Path_CRP_Curves_ShankFoot = [Output_File_Path 'CRP_Curve_ShankFoot' '.xls']; 

xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Knee_Two, 'Knee 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Knee_Four, 'Knee 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Knee_Six, 'Knee 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Knee_Eight, 'Knee 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Hip_Two, 'Hip 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Hip_Four, 'Hip 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Hip_Six, 'Hip 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Hip_Eight, 'Hip 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Umbilical_Two, 'Umbilical 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Umbilical_Four, 'Umbilical 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Umbilical_Six, 'Umbilical 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Umbilical_Eight, 'Umbilical 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Two, 'Xiphoid 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Four, 'Xiphoid 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Six, 'Xiphoid 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_ShankFoot, CRP_Curves_ShankFoot_Export_Mean_Xiphoid_Eight, 'Xiphoid 0.8'); 


% Export CRP Curve 0-100% ThighShank
Output_File_Path_CRP_Curves_ThighShank = [Output_File_Path 'CRP_Curve_ThighShank' '.xls']; 

xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Knee_Two, 'Knee 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Knee_Four, 'Knee 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Knee_Six, 'Knee 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Knee_Eight, 'Knee 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Hip_Two, 'Hip 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Hip_Four, 'Hip 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Hip_Six, 'Hip 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Hip_Eight, 'Hip 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Umbilical_Two, 'Umbilical 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Umbilical_Four, 'Umbilical 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Umbilical_Six, 'Umbilical 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Umbilical_Eight, 'Umbilical 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Xiphoid_Two, 'Xiphoid 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Xiphoid_Four, 'Xiphoid 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Xiphoid_Six, 'Xiphoid 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_ThighShank, CRP_Curves_ThighShank_Export_Mean_Xiphoid_Eight, 'Xiphoid 0.8'); 



% Export CRP Curve 0-100% TrunkThigh
Output_File_Path_CRP_Curves_TrunkThigh = [Output_File_Path 'CRP_Curve_TrunkThigh' '.xls']; 

xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Knee_Two, 'Knee 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Knee_Four, 'Knee 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Knee_Six, 'Knee 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Knee_Eight, 'Knee 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Hip_Two, 'Hip 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Hip_Four, 'Hip 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Hip_Six, 'Hip 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Hip_Eight, 'Hip 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Two, 'Umbilical 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Four, 'Umbilical 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Six, 'Umbilical 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Umbilical_Eight, 'Umbilical 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Two, 'Xiphoid 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Four, 'Xiphoid 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Six, 'Xiphoid 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_TrunkThigh, CRP_Curves_TrunkThigh_Export_Mean_Xiphoid_Eight, 'Xiphoid 0.8'); 

% Export CRP Curve 0-100% KneeAnkle
Output_File_Path_CRP_Curves_KneeAnkle = [Output_File_Path 'CRP_Curve_KneeAnkle' '.xls']; 

xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Knee_Two, 'Knee 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Knee_Four, 'Knee 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Knee_Six, 'Knee 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Knee_Eight, 'Knee 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Hip_Two, 'Hip 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Hip_Four, 'Hip 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Hip_Six, 'Hip 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Hip_Eight, 'Hip 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Two, 'Umbilical 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Four, 'Umbilical 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Six, 'Umbilical 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Umbilical_Eight, 'Umbilical 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Two, 'Xiphoid 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Four, 'Xiphoid 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Six, 'Xiphoid 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_KneeAnkle, CRP_Curves_KneeAnkle_Export_Mean_Xiphoid_Eight, 'Xiphoid 0.8'); 

%HipKnee
Output_File_Path_CRP_Curves_HipKnee = [Output_File_Path 'CRP_Curve_HipKnee' '.xls']; 

xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Knee_Two, 'Knee 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Knee_Four, 'Knee 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Knee_Six, 'Knee 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Knee_Eight, 'Knee 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Hip_Two, 'Hip 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Hip_Four, 'Hip 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Hip_Six, 'Hip 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Hip_Eight, 'Hip 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Umbilical_Two, 'Umbilical 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Umbilical_Four, 'Umbilical 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Umbilical_Six, 'Umbilical 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Umbilical_Eight, 'Umbilical 0.8'); 

xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Xiphoid_Two, 'Xiphoid 0.2'); 
xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Xiphoid_Four, 'Xiphoid 0.4'); 
xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Xiphoid_Six, 'Xiphoid 0.6'); 
xlswrite(Output_File_Path_CRP_Curves_HipKnee, CRP_Curves_HipKnee_Export_Mean_Xiphoid_Eight, 'Xiphoid 0.8'); 





%% Figures

% figure ('Name','Segments Angles')
% title ('Segments Angles');
% xlabel('Time (s)');
% ylabel ('Angle (°)');
% hold on
% plot (Foot_Angle_Stride_1, 'r', 'LineWidth', 2);
% hold on
% plot (Shank_Angle_Stride_1, 'b', 'LineWidth', 2);
% hold on
% legend ('Foot', 'Shank')
% 
% figure ('Name','Segments Angles Centered')
% title ('Segments Angles Centered');
% xlabel('Time (s)');
% ylabel ('Angle (°)');
% hold on
% plot (Foot_Angle_Stride_1_Centered, 'r', 'LineWidth', 2);
% hold on
% plot (Shank_Angle_Stride_1_Centered, 'b', 'LineWidth', 2);
% hold on
% legend ('Foot', 'Shank')
% 
% figure ('Name','Segments Angles Hilbert')
% title ('Segments Angles Hilbert');
% xlabel('Time (s)');
% ylabel ('Angle (°)');
% hold on
% plot (Foot_Angle_Stride_1_Hilbert, 'r', 'LineWidth', 2);
% hold on
% plot (Shank_Angle_Stride_1_Hilbert, 'b', 'LineWidth', 2);
% hold on
% legend ('Foot', 'Shank')
% 
% 
% figure ('Name','Segments Angles Hilbert')
% title ('Phase Angle');
% xlabel('Time (s)');
% ylabel ('Angle (°)');
% hold on
% plot (Phase_Angle_Foot_Stride_1, 'r', 'LineWidth', 2);
% hold on
% plot (Phase_Angle_Shank_Stride_1, 'b', 'LineWidth', 2);
% hold on
% legend ('Foot', 'Shank')
% 
% figure ('Name','CRP Shank Foot')
% title ('CRP Shank Foot');
% xlabel('Time (s)');
% ylabel ('Angle (°)');
% hold on
% plot (CRP_ShankFoot_Stride_1, 'r', 'LineWidth', 2);
% hold on
% legend ('CRP Shank Foot')
