%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef Args < handle
    properties(SetAccess=private)
        p;
        argued;
    end
    
    methods
        function this=Args(p, varargin)
            this.p=p;
            parse(p,varargin{:});
            this.argued=SortedStringSet.New(p);
        end
    end
    methods(Static)
        function [args, argued, unmatchedArgs]=NewKeepUnmatched(p, varargin)
            p.KeepUnmatched=true;
            this=Args(p, varargin{:});
            args=this.p.Results;
            argued=this.argued;
            unmatchedArgs=this.p.Unmatched;
            if isempty(fieldnames(unmatchedArgs))
                unmatchedArgs=[];
            end
        end

        function [args, argued]=New(p, varargin)
            this=Args(p, varargin{:});
            args=this.p.Results;
            argued=this.argued;
        end
        
        function argument=Get(name, varargin)
            argument=[];
            N=length(varargin);
            for i=1:2:N
                if strcmpi(name, varargin{i})
                    argument=varargin{i+1};
                    break;
                end
            end
        end
        
        function varargin=Set(name, value, varargin)
            N=length(varargin);
            for i=1:2:N
                if strcmpi(name, varargin{i})
                    varargin{i+1}=value;
                    return;
                end
            end
            varargin{end+1}=name;
            varargin{end+1}=value;
        end
    end
end