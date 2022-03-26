        subroutine upcase(string)
                                 
c       Translates a string to upper case.
                                          
        character*(*) string
                            
        integer i
                 
        do 310 i=1, len(string)
          if (ichar(string(i:i)).ge.97 .and. ichar(string(i:i)).le.122)
     &    string(i:i) = char(ichar(string(i:i))-32)
310     continue
                
        end
