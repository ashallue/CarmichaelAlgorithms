


bool Preproduct::is_admissable( int64_t prime_to_append )
{
  while( prime_to_append > next_inadmissable[0] )
  {
    next_inadmissable[0] += (appended_primes[0] << mod_three_status[0]);
    mod_three_status[0] = (mod_three_status[0] == 1 ) ? 2 : 1;
    int i= 1;
    while( i < len_appended_primes && next_inadmissable[ i-1 ] < next_inadmissable[ i ] )
    {
      std::swap( next_inadmissable[ i-1 ], next_inadmissable[ i ] );
      std::swap( appended_primes[ i-1 ], appended_primes[ i ] );
      std::swap( mod_three_status[ i-1 ], mod_three_status[ i ] );
      i++;
    }
  }
  return ( prime_to_append < next_inadmissable[0] ) ;
}
