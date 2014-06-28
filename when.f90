        ! By Haiyun Jin 2014-1-15
        ! return current system time
      
      subroutine when
        implicit none

      integer*4 today(3), now(3)

      call idate(today)   ! today(1)=day, (2)=month, (3)=year
      call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
      write ( *, 100 )  today(2), today(1), today(3), now
 100  format ( 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ', i2.2, ':', i2.2, ':', i2.2 )

       end  subroutine   ! when
