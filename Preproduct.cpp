


bool Preproduct::is_admissible( int64_t prime_to_append )
{
  // if this works how we want it to, this while loop will not be entered often
  while( prime_to_append > next_inadmissible[0] )
  {
    // add 2*p or 4*p and avoid divisibility be 3.  Flip the state of mod 3 status.
    next_inadmissible[0] += (appended_primes[0] << mod_three_status[0]);
    mod_three_status[0] = (mod_three_status[0] == 1 ) ? 2 : 1;
    int i= 1;
    // put next_inadmissible in sorted order having increased the first element
    while( i < len_appended_primes && next_inadmissible[ i-1 ] < next_inadmissible[ i ] )
    {
      std::swap( next_inadmissible[ i-1 ], next_inadmissible[ i ] );
      std::swap( appended_primes[ i-1 ], appended_primes[ i ] );
      std::swap( mod_three_status[ i-1 ], mod_three_status[ i ] );
      i++;
    }
  }
  // now prime to append is less than or equal to next_inadmissible
  // return true if it is less than; return false if it is equal
  return ( prime_to_append < next_inadmissible[0] ) ;
}
