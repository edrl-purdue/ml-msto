classdef parfor_waitbar < handle
    %This class creates a waitbar or message when using for or parfor.
    %Required Input:
    %TotalMessage: N in "i = 1: N".
    %Optional Inputs:
    %'Waitbar': true of false (default). If true, this class creates a
    %               waitbar.
    %'FileName': 'screen' or a char array. If 'screen', print the message
    %               on screen; otherwise, save the message in the file
    %               named 'FileName'.
    %'ReportInterval': 1x1. Report at every i is costly. This number
    %                defines the interval for reporting. 
    %%To use this class, one needs to call the class right before the loop:
    %N = 1000;
    %WaitMessage = parfor_waitbar(N);
    %%Call "Send" method in the loop. 
    %for i = 1: N
    %   WaitMessage.Send;
    %   pause(0.5);
    %end
    %%Delete the obj after the loop.
    %WaitMessage.Destroy;
    %Copyright (c) 2019, Yun Pu
    properties (SetAccess = private)
        NumMessage; %Number of messages received from the workers.
        TotalMessage; %Number of total messages.
        Waitbar; %If waitbar is true, create a waitbar; otherwise, save the message in a file.
        FileName; %If FileName = 'screen', the current message does not save in a file.
        StartTime
        UsedTime_1; %Time at last step.
        WaitbarHandle;
        ReportInterval;
        FileID;
        DataQueueHandle;
        bar0;
        bar10;
        bar20;
        bar30;
        bar40;
        bar50;
        bar60;
        bar70;
        bar80;
        bar90;
        bar100;
    end
    
   methods
       function Obj = parfor_waitbar(TotalMessage, varargin)
           Obj.DataQueueHandle = parallel.pool.DataQueue;
           Obj.StartTime = tic;
           Obj.NumMessage = 0;
           Obj.UsedTime_1 = Obj.StartTime;
           Obj.TotalMessage = TotalMessage;
           Obj.bar0 = 0;
           Obj.bar10 = 0;
           Obj.bar20 = 0;
           Obj.bar30 = 0;
           Obj.bar40 = 0;
           Obj.bar50 = 0;
           Obj.bar60 = 0;
           Obj.bar70 = 0;
           Obj.bar80 = 0;
           Obj.bar90 = 0;
           Obj.bar100 = 0;
           InParser = inputParser;
           addParameter(InParser,'Waitbar',false,@islogical);
           addParameter(InParser,'FileName', 'screen', @ischar);
           addParameter(InParser,'ReportInterval', ceil(TotalMessage/100), @isnumeric);
           parse(InParser, varargin{:})
           Obj.Waitbar = InParser.Results.Waitbar;
           Obj.FileName = InParser.Results.FileName;
           Obj.ReportInterval = InParser.Results.ReportInterval;
           if Obj.Waitbar
               Obj.WaitbarHandle = waitbar(0, [num2str(0), '%'], 'Resize', true);
           end
           switch Obj.FileName
               case 'screen'
               otherwise
                   Obj.FileID = fopen(Obj.FileName, 'w');     
           end       
           afterEach(Obj.DataQueueHandle, @Obj.Update);
       end
       function Send(Obj)
           send(Obj.DataQueueHandle, 0);
       end
       function Destroy(Obj)
           if Obj.Waitbar
               delete(Obj.WaitbarHandle);
           end
           delete(Obj.DataQueueHandle);
           delete(Obj);
       end
   end
   
   methods (Access = private)
       function Obj = Update(Obj, ~)
           Obj.AddOne;
           if mod(Obj.NumMessage, Obj.ReportInterval)
               return
           end
           if Obj.Waitbar             
               Obj.WaitbarUpdate;
           else
               Obj.FileUpdate;
           end
       end
       
       function WaitbarUpdate(Obj)
           UsedTime_now = toc(Obj.StartTime);
           EstimatedTimeNeeded = (UsedTime_now-Obj.UsedTime_1)/Obj.ReportInterval*(Obj.TotalMessage-Obj.NumMessage);
           waitbar(Obj.NumMessage/Obj.TotalMessage, Obj.WaitbarHandle, [num2str(Obj.NumMessage/Obj.TotalMessage*100, '%.2f'), '%; ', num2str(UsedTime_now, '%.2f'), 's used and ', num2str(EstimatedTimeNeeded, '%.2f'), 's needed.']);
           Obj.UsedTime_1 = UsedTime_now;
       end
       
       function FileUpdate(Obj)
           UsedTime_now = toc(Obj.StartTime);
           EstimatedTimeNeeded = (UsedTime_now-Obj.UsedTime_1)/Obj.ReportInterval*(Obj.TotalMessage-Obj.NumMessage);           
           switch Obj.FileName
               case 'screen'
%                    fprintf('%.2f%%; %.2fs used and %.2fs needed...\n', Obj.NumMessage/Obj.TotalMessage*100, UsedTime_now, EstimatedTimeNeeded);
                   if Obj.NumMessage/Obj.TotalMessage*100 > 0 && Obj.bar0 == 0
                       fprintf('     |=');
                       Obj.bar0 = 1;
                   end
                   if Obj.NumMessage/Obj.TotalMessage*100 > 10 && Obj.bar10 == 0
                       fprintf('=====');
                       Obj.bar10 = 1;
                   end
                   if Obj.NumMessage/Obj.TotalMessage*100 > 20 && Obj.bar20 == 0
                       fprintf('=====');
                       Obj.bar20 = 1;
                   end
                   if Obj.NumMessage/Obj.TotalMessage*100 > 30 && Obj.bar30 == 0
                       fprintf('=====');
                       Obj.bar30 = 1;
                   end
                   if Obj.NumMessage/Obj.TotalMessage*100 > 40 && Obj.bar40 == 0
                       fprintf('=====');
                       Obj.bar40 = 1;
                   end
                   if Obj.NumMessage/Obj.TotalMessage*100 > 50 && Obj.bar50 == 0
                       fprintf('=====');
                       Obj.bar50 = 1;
                   end
                   if Obj.NumMessage/Obj.TotalMessage*100 > 60 && Obj.bar60 == 0
                       fprintf('=====');
                       Obj.bar60 = 1;
                   end
                   if Obj.NumMessage/Obj.TotalMessage*100 > 70 && Obj.bar70 == 0
                       fprintf('=====');
                       Obj.bar70 = 1;
                   end
                   if Obj.NumMessage/Obj.TotalMessage*100 > 80 && Obj.bar80 == 0
                       fprintf('=====');
                       Obj.bar80 = 1;
                   end
                   if Obj.NumMessage/Obj.TotalMessage*100 > 90 && Obj.bar90 == 0
                       fprintf('=====');
                       Obj.bar90 = 1;
                   end
                   if Obj.NumMessage/Obj.TotalMessage*100 == 100 && Obj.bar100 == 0
                       fprintf('=====|\n');
                       Obj.bar100 = 1;
                   end
               otherwise
                   fprintf(Obj.FileID, '%.2f%%; %.2fs used and %.2fs needed...\n', Obj.NumMessage/Obj.TotalMessage*100, UsedTime_now, EstimatedTimeNeeded);
           end
           Obj.UsedTime_1 = UsedTime_now;
       end
       function AddOne(Obj)
           Obj.NumMessage = Obj.NumMessage + 1;
       end      
   end
end