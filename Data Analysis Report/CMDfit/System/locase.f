        subroutine locase(string)
                                 
c       Translates a string to lower case.
                                          
        character*(*) string
                            
        integer i
                 
        do 310 i=1, len(string)
          if (ichar(string(i:i)).ge.65 .and. ichar(string(i:i)).le.90)
     &    string(i:i) = char(ichar(string(i:i))+32)
310     continue
                
        end
