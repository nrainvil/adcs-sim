function ym = monthofyear( varargin )
%MONTHOFYEAR Ordinal number of month in year.
%
%   MONTHOFYEAR(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND) returns the ordinal
%   month number in the given year plus a fractional part depending on the
%   day and time of day.
%
%   Any missing MONTH or DAY will be replaced by ones.  Any missing HOUR,
%   MINUTE or SECOND will be replaced by zeros.
%
%   If no date is specified, the current date and time is used.  Gregorian
%   calendar is assumed.

%   Author:      Peter John Acklam
%   Time-stamp:  2002-03-03 12:50:16 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   nargsin = nargin;
   error(nargchk(0, 6, nargsin));
   if nargsin
      argv = {1 1 1 0 0 0};
      argv(1:nargsin) = varargin;
   else
      argv = num2cell(clock);
   end
   [year, month, day, hour, minute, second] = deal(argv{:});

   ym = month + (dayofmonth(year, month, day, hour, minute, second)) ...
                ./ daysinmonth(year, month);