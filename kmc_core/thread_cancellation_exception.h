#pragma once

//In case of critical error in one thread
//the information of critical situation is send to other threads
//which are throwing the exception of below class
//so it is just a information that thread is finishing because of critical situation in other thread
//the CExceptionAwareThread will ignore this exception, but the real std::thread will end
class CThreadCancellationException { };
