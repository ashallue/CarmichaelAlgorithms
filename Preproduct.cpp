


bool Preproduct::is_admissable( int64_t prime_to_append )
{
  bool return_val;
  if( primes_appended == 0 || prime_to_append < next_inadmissable[0] )
  {
    return_val = true;
  }
  else // primes_appended > 0 and prime_to_append >= next_inadmissible[0]
  {
    if( prime_to_append == next_inadmissable[0] )
    {
        return_val = false;
    }
    else // primes_appended > 0 and prime_to_append > next_inadmissible[0]
    {
      // this could be implemented in a while loop
      // however, the recursive call will force that to happen
      // we take the risk that a single iteration will usually suffice
      next_inadmissible[0] += (appended_primes[0] << (mod_three_status[0] + 1) )
      mod_three_status[0] = ( mod_three_status[0] + 1 ) % 2;

      if( primes_appended > 1 )
      {
        
        // need to update
        // next_inadmissible_prime[5];
        // bool mod_three_status[5];
        // uint64_t P_primes[5];
        // via bubble sort
      }

      return_val =  Preproduct::is_admissable( prime_to_append );

    }
  }
  return return_val;
}
