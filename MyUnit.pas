Unit myunit; 
{--------------------------------} 
Interface 
Type 
	Complex = record 
		Re, Im: real; 
	end; 
	TMassExtended = array of extended;

Procedure AddC( x, y: complex; var z: complex );
Procedure ReadMVF( str: string; var m: integer; var Mass, Vector: TMassExtended );
Procedure PrintMatrix( m: integer; Mass: TMassExtended );
Procedure PrintVector( m: integer; Vector: TMassExtended );
Procedure MultMV( m: integer; Matx, Vector: TMassExtended; var v0: TMassExtended );
Procedure DivisionVV( m: integer; v1, v2: TMassExtended; var v0: TMassExtended );
Procedure Normalize( m: integer; var v: TMassExtended );
Procedure PrintMaxElVector( m: integer; v: TMassExtended );
Function FinalCheck ( m: integer; v1, v0: TMassExtended ):boolean;
Function Norma( m: integer; v: TMassExtended ):extended;
Procedure SwapVectorMass( m: integer; var A: TMassExtended; k, n: integer );
Procedure MethodGauss( m: integer; A: TMassExtended; b: TMassExtended; var x: TMassExtended; var erorr: byte );
{---------------------------------} 
Implementation 
	Procedure AddC( x, y: complex; var z: complex ); 
	begin
		z.re:= x.re + y.re; 
		z.im:= x.im + y.im; 
	end ; 
	
	Procedure ReadMVF( str: string; var m: integer; var Mass, Vector: TMassExtended );
	var f: textfile;
		i: integer;
	begin
		assign( f, str );
		reset( f );
			read( f, m );
			setlength( Mass, m * m );
			setlength( Vector, m );
			for i:= 0 to m * m - 1 do
				read( f, Mass[ i ] );
			for i:= 0 to m - 1 do
				read( f, Vector[ i ] );
		close( f );
	end;
	
	Procedure PrintMatrix( m: integer; Mass: TMassExtended );
	var	i, j: integer;
	begin
		for i:= 0 to m - 1 do
			begin
				for j:= 0 to m - 1 do
					write( Mass[ i * m + j ]:8:8, ' ' );
				writeln;
			end;
	end;
	
	Procedure PrintVector( m: integer; Vector: TMassExtended );
	var i: integer;
	begin
		for i:= 0 to m - 1 do
			write( Vector[ i ], ' ' );
		writeln;
	end;
	
	Procedure MultMV( m: integer; Matx, Vector: TMassExtended; var v0: TMassExtended );
	var i, j: integer;
		sum: extended;
	begin
		setlength( v0, m );
		for i:= 0 to m - 1 do
			begin
				sum:= 0;
				for j:= 0 to m - 1 do	
					sum:= sum + Matx[ i * m + j ] * Vector[ j ];
				v0[ i ]:= sum;
			end;
				
	end;
	
	Procedure DivisionVV( m: integer; v1, v2: TMassExtended; var v0: TMassExtended );
	var i: integer;
		eps: real;
	begin
		eps:= 0.000001;
		setlength( v0, m);
		for i:= 0 to m - 1 do
			if abs( v2[ i ] ) > eps then v0[ i ]:= v1[ i ] / v2[ i ];
	end;
	
	Function FinalCheck( m: integer; v1, v0: TMassExtended ):boolean;
	var	i: integer;
		bool: boolean;
		eps: real;
	begin
		eps:= 0.000000001;
		for i:= 0 to m - 1 do
			if abs( v1[ i ] - v0[ i ] ) < eps then  begin
				if i = m - 1 then begin 
					bool:= true;
					break;
				end;
			end
			else bool:= false;
		
		FinalCheck:= bool;
	end;
	
	Function Norma( m: integer; v: TMassExtended ):extended;
	var i: integer;	
		sum: extended;
	begin
		sum:= 0;
		for i:= 0 to m - 1 do
			sum:= sum + v[ i ] * v[ i ];
		Norma:= sqrt ( sum );
	end;
	
	Procedure Normalize( m: integer; var v: TMassExtended );
	var i: integer;
		ex, eps: extended;
	begin
		setlength( v, m );
		ex:= norma( m, v );
		eps:= 0.000001;
		if ex > eps then
			for i:= 0 to m - 1 do
				v[ i ]:= v[ i ] / ex;
	end;
	
	Procedure PrintMaxElVector( m: integer; v: TMassExtended );
	var i, k: integer;
	begin
		k:= 0;
		for i:= 0 to m - 1 do
			if abs( v[ i ] ) > abs( v[ k ] ) then k:= i;
		writeln( v[ k ] );
	end;
	
	Procedure SwapVectorMass( m: integer; var A: TMassExtended; k, n: integer );
	var ext: extended;
		i: integer;
	begin
		for i:= 0 to m - 1 do
			begin
				ext:= A[ k * m + i ];
				A[ k * m + i ]:= A[ n * m + i ];
				A[ n * m + i ]:= ext;
			end;
	end;
	
	Procedure MethodGauss( m: integer; A: TMassExtended; b: TMassExtended; var x: TMassExtended; var erorr: byte );
	var i, j, k: integer;
		ext, c: extended;
	begin
		setlength( x, m );
		for i:= 0 to m - 2 do
			begin
				if A[ i * m + i ] = 0 then 
					begin
						for j:= i + 1 to m - 1 do
							begin
								if A[ j * m + i ] <> 0 then
									SwapVectorMass( m, A, j, i)
								else 
									begin
										if j = ( m - 1 ) then 
											erorr:= 100;
									end;
							end;
						if erorr = 100 then break;
					end;
					
				for j:= i + 1 to m - 1 do 
					begin
						c:= A[ j * m + i ] / A[ i * m + i ];
						for k:= i to m - 1 do
							begin
								A[ j * m + k ]:= A[ j * m + k ] - A[ i * m + k ] * c;
							end;
						b[ j ]:= b[ j ] - b[ i ] * c; 
					end;
			end;
			printmatrix(m, A);
			printvector( m, b);
		// Обратный ход
		for i:= ( m - 1 ) downto 0 do
			begin
				if A[ i * m + i ] = 0 then begin
					erorr:= 101;
					break;
				end
				else
					begin
						c:= 0;
						for j:= i + 1 to m - 1 do 
							c:= c + x[ j ] * A[ i * m + j ];
							write( c, ' ' );
						x[ i ]:= ( b[ i ] - c ) / A[ i * m + i ];
						writeln;
					end;
			end;
	end;
end.