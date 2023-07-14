%% 
% ================= ANÁLISE SACCADAS =================

%%
clear all; close all; clc;
addpath(fullfile('..','tools','edfImport'))

%% 
group = "P";    

%% 
for s = 1
    
    if s < 10
        sub = strcat(group,'0',num2str(s));
    else
        sub = strcat(group,num2str(s));
    end
    
    dir_files = strcat('F:\SEM_mri\MRI\SEM_mri_rawdata\', sub, '\ET\'); % Directoria
    
    runs = [1:7];
    
    for r = 1%runs
        
        nrun = strcat('R',num2str(r));
        edfFile = strcat(sub, '_', nrun, '.edf');
        
        %% PARÂMETROS
        
        % Dimensões ecrã
        screenXpixels = 1920; 
        screenYpixels = 1080; 
        
        eye = 1; % Olho com que o eye tracker foi adquirido
        
        % Limites horizontais a partir dos quais se deteta uma sacada para
        % a direita ou esquerda (relacionados com as linhas verticais à
        % esquerda e direita presentes no paradigma)
        left_line = 200; 
        right_line = screenXpixels - 200;
        dist_lim = 320;
        left_lim = left_line+dist_lim;
        right_lim = right_line-dist_lim;
                
        % Distância aos limites verticais do ecrã até aos quais se pode
        % detetar uma resposta (Ys elevados poderão ser blinks/falhas, por ex)
        dist_ylim = 200;
        
        % Sampling Rate
        srate = 1000;
                
        % Nº mínimo de pontos para considerar que houve movimento ocular
        % para a direita ou esquerda
        points_min = 15;
        
        % Valor de x a partir do qual se considera a deteção de movimento ocular
        xmin = 0;

        %% Data

        run = edfImport(char(strcat(dir_files,"\" ,edfFile)), [1 1 1], '');
        triggers = load(strcat('F:\SEM_mri\MRI\SEM_mri_rawdata\', sub, '\mat_files\', ...
            sub, '_', nrun, '.mat'));
        
        %% Triggers

        instTime = triggers.Output.Run.dataMat(6,:);
        responseTime = triggers.Output.Run.dataMat(8,:);
        gap2Time = triggers.Output.Run.dataMat(9,:);
        noiseTime = triggers.Output.Run.NoiseTime;

        %% Análise das respostas

        Response = {};
        trials = [1:length(triggers.Output.Run.dataMat)]; 
        close all
        figure
        for i = trials 

            %==== Coordenadas oculares ====

            % Gaze during response time
            xGazeResponse = run.Samples.gx(eye,round(responseTime(i),1)*srate:round(gap2Time(i),1)*srate-1);
            yGazeResponse = run.Samples.gy(eye,round(responseTime(i),1)*srate:round(gap2Time(i),1)*srate-1);

            % Gaze before response time
            xGazeBefore = run.Samples.gx(eye,round(instTime(i)+0.05,1)*srate:round(responseTime(i),1)*srate-1);
            yGazeBefore = run.Samples.gy(eye,round(instTime(i)+0.05,1)*srate:round(responseTime(i),1)*srate-1);

            % Gaze after response time
            if i < length(triggers.Output.Run.dataMat) && ...
                    sum(ismember(round(noiseTime),[round(gap2Time(i))+1, round(gap2Time(i))+2, round(gap2Time(i))+3])) == 0
                xGazeAfter = run.Samples.gx(eye,round(gap2Time(i),1)*srate:round(instTime(i+1),1)*srate-1);
                yGazeAfter = run.Samples.gy(eye,round(gap2Time(i),1)*srate:round(instTime(i+1),1)*srate-1);
            elseif i == length(triggers.Output.Run.dataMat)
                xGazeAfter = run.Samples.gx(eye,round(gap2Time(i),1)*srate:round(noiseTime(3),1)*srate-1);
                yGazeAfter = run.Samples.gy(eye,round(gap2Time(i),1)*srate:round(noiseTime(3),1)*srate-1);
            elseif isequal(ismember(round(noiseTime), [round(gap2Time(i))+1, round(gap2Time(i))+2, round(gap2Time(i))+3]), [1 0 0])
                xGazeAfter = run.Samples.gx(eye,round(gap2Time(i),1)*srate:round(noiseTime(1),1)*srate-1);
                yGazeAfter = run.Samples.gy(eye,round(gap2Time(i),1)*srate:round(noiseTime(1),1)*srate-1);
            else
                xGazeAfter = run.Samples.gx(eye,round(gap2Time(i),1)*srate:round(noiseTime(2),1)*srate-1);
                yGazeAfter = run.Samples.gy(eye,round(gap2Time(i),1)*srate:round(noiseTime(2),1)*srate-1);
            end
            
            dist_pix = 500;
            
            outGazeResponse = find(xGazeResponse < -dist_pix | xGazeResponse > screenXpixels + dist_pix | ...
                yGazeResponse < -dist_pix | yGazeResponse > screenYpixels + dist_pix); 
            xGazeResponse(outGazeResponse) = [];
            yGazeResponse(outGazeResponse) = [];

            outGazeBefore = find(xGazeBefore < -dist_pix | xGazeBefore > screenXpixels + dist_pix | ...
                yGazeBefore < -dist_pix | yGazeBefore > screenYpixels + dist_pix); 
            xGazeBefore(outGazeBefore) = [];
            yGazeBefore(outGazeBefore) = [];

            outGazeAfter = find(xGazeAfter < -dist_pix | xGazeAfter > screenXpixels + dist_pix | ...
                yGazeAfter < -dist_pix | yGazeAfter > screenYpixels + dist_pix); 
            xGazeAfter(outGazeAfter) = [];
            yGazeAfter(outGazeAfter) = [];

            %==== Cálculo do sentido da sacada ====

            % Tipos de resposta:
            % Direita (antes, durante e/ou depois);
            % Esquerda (antes, durante e/ou depois);
            % Ambos (antes, durante e/ou depois);
            % Fixação.

            % Pre-response time
            if length(find(xGazeBefore>right_lim & yGazeBefore>dist_ylim & yGazeBefore <screenYpixels-dist_ylim)) > points_min && ...
                    length(find(xGazeBefore>xmin & xGazeBefore<left_lim & yGazeBefore>dist_ylim & yGazeBefore <screenYpixels-dist_ylim)) > points_min
                idx_right = find(xGazeBefore>right_lim & yGazeBefore>dist_ylim & yGazeBefore <screenYpixels-dist_ylim);
                idx_left = find(xGazeBefore>xmin & xGazeBefore<left_lim & yGazeBefore>dist_ylim & yGazeBefore <screenYpixels-dist_ylim);
                if idx_right(1)==1 && sum(diff(idx_right)~=1) == 0 
                    Response{i,1} = {'Left'};
                elseif idx_left(1)==1 && sum(diff(idx_left)~=1) == 0 
                    Response{i,1} = {'Right'};
                else
                    Response{i,1} = {'Both'};
                end
            elseif length(find(xGazeBefore>right_lim & yGazeBefore>dist_ylim & yGazeBefore <screenYpixels-dist_ylim)) > points_min
                % Verificar 1º se não é uma resposta tardia à trial anterior ("volta" com o sentido contrário)
                idx = find(xGazeBefore>right_lim & yGazeBefore>dist_ylim & yGazeBefore <screenYpixels-dist_ylim);
                if idx(1)==1 && sum(diff(idx)~=1) == 0 
                    %Se o 1º ponto do dataset já estiver à direita e se não houver um regresso ao lado direito, significa que é apenas uma resposta tardia à trial anterior
                    Response{i,1} = {'Fixation'};
                else
                    Response{i,1} = {'Right'};
                end
            elseif length(find(xGazeBefore>xmin & xGazeBefore<left_lim & yGazeBefore>dist_ylim & yGazeBefore <screenYpixels-dist_ylim)) > points_min
                % Verificar 1º se não é uma resposta tardia à trial anterior ("volta" com o sentido contrário)
                idx = find(xGazeBefore>xmin & xGazeBefore<left_lim & yGazeBefore>dist_ylim & yGazeBefore <screenYpixels-dist_ylim);
                if idx(1)==1 && sum(diff(idx)~=1) == 0 
                    Response{i,1} = {'Fixation'};
                else
                    Response{i,1} = {'Left'};
                end
%             elseif length(xGazeBefore) ==0
%                 Response{i,1} = {'None'};
            else
                Response{i,1} = {'Fixation'};
            end

            % Response time
            if length(find(xGazeResponse>right_lim & yGazeResponse>dist_ylim & yGazeResponse <screenYpixels-dist_ylim)) > points_min && ...
                    length(find(xGazeResponse>xmin & xGazeResponse<left_lim & yGazeResponse>dist_ylim & yGazeResponse <screenYpixels-dist_ylim)) > points_min
                Response{i,2} = {'Both'};
            elseif length(find(xGazeResponse>right_lim & yGazeResponse>dist_ylim & yGazeResponse <screenYpixels-dist_ylim)) > points_min
                Response{i,2} = {'Right'};
            elseif length(find(xGazeResponse>xmin & xGazeResponse<left_lim & yGazeResponse>dist_ylim & yGazeResponse <screenYpixels-dist_ylim)) > points_min
                Response{i,2} = {'Left'};
%             elseif length(xGazeResponse) ==0
%                 Response{i,1} = {'None'};
            else
                Response{i,2} = {'Fixation'};
            end

            % Post-response time
            if length(find(xGazeAfter>right_lim & yGazeAfter>dist_ylim & yGazeAfter <screenYpixels-dist_ylim)) > points_min && ...
                    length(find(xGazeAfter>xmin & xGazeAfter<left_lim & yGazeAfter>dist_ylim & yGazeAfter <screenYpixels-dist_ylim)) > points_min
                Response{i,3} = {'Both'};
            elseif length(find(xGazeAfter>right_lim & yGazeAfter>dist_ylim & yGazeAfter <screenYpixels-dist_ylim)) > points_min
                Response{i,3} = {'Right'};
            elseif length(find(xGazeAfter>xmin & xGazeAfter<left_lim & yGazeAfter>dist_ylim & yGazeAfter <screenYpixels-dist_ylim)) > points_min
                Response{i,3} = {'Left'};
%             elseif length(xGazeAfter) ==0
%                 Response{i,1} = {'None'};
            else
                Response{i,3} = {'Fixation'};
            end

            %% 

            %==== Figuras ====

            xmin_plot = 0;
            xmax_plot = screenXpixels;
            ymin_plot = 0;
            ymax_plot = screenYpixels;

            plot(xGazeResponse, screenYpixels-yGazeResponse,'.r','LineWidth',2)
            xline(left_lim, '--')
            xline(right_lim, '--')
            xlim([xmin_plot xmax_plot])
            ylim([ymin_plot ymax_plot])
            title('Response', 'fontsize', 18)
            if length(trials)>1
                hold on
            else
                figure
                plot(xGazeBefore, screenYpixels-yGazeBefore,'.r','LineWidth',2)
                xline(left_lim, '--')
                xline(right_lim, '--')
                xlim([xmin_plot xmax_plot])
                ylim([ymin_plot ymax_plot])
                title('Pre-response', 'fontsize', 18)

                figure
                plot(xGazeAfter, screenYpixels-yGazeAfter,'.r','LineWidth',2)
                xline(left_lim, '--')
                xline(right_lim, '--')
                xlim([xmin_plot xmax_plot])
                ylim([ymin_plot ymax_plot])
                title('Post-response', 'fontsize', 18)
            end

            xlim([xmin_plot xmax_plot])
            ylim([ymin_plot ymax_plot])

        end
    end
end
