#pragma once

#include "data_structure/parallel/cache.h"
#include "data_structure/parallel/metaprogramming_utils.h"
#include "data_structure/parallel/spin_lock.h"

#include <algorithm>
#include <atomic>
#include <functional>
#include <future>
#include <memory>
#include <thread>
#include <type_traits>
#include <vector>

#ifdef __gnu_linux__
#include <numa.h>
#include <pthread.h>
#endif

namespace parallel {

static void PinToCore(size_t core) {
#ifdef __gnu_linux__
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core, &cpuset);
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#endif
}

static void Unpin() {
#ifdef __gnu_linux__
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    for (size_t core = 0; core < std::thread::hardware_concurrency(); ++core) {
        CPU_SET(core, &cpuset);
    }
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#endif
}

template<typename T>
class TThreadsafeQueue {
private:
        struct TNode {
                std::shared_ptr <T> Data;
                std::unique_ptr <TNode> Next;
        };

        parallel::spin_lock HeadMutex;
        std::unique_ptr <TNode> Head;
        parallel::spin_lock TailMutex;
        TNode* Tail;

        TNode* GetTail() {
                std::lock_guard <parallel::spin_lock> tailLock(TailMutex);
                return Tail;
        }

        bool PopHead(T& value) {
                std::lock_guard <parallel::spin_lock> headLock(HeadMutex);

                if (Head.get() == GetTail()) {
                        return false;
                }

                value = std::move(*Head->Data);
                std::unique_ptr <TNode> const oldHead = std::move(Head);
                Head = std::move(oldHead->Next);
                return oldHead.get() != nullptr;
        }

        template<typename TType>
        void PushImpl(TType&& newValue) {
                std::shared_ptr <T> newData(std::make_shared<T>(std::forward<TType>(newValue)));
                std::unique_ptr <TNode> p(new TNode);
                TNode* const newTail = p.get();
                std::lock_guard <parallel::spin_lock> tailLock(TailMutex);
                Tail->Data = newData;
                Tail->Next = std::move(p);
                Tail = newTail;
        }


public:
        TThreadsafeQueue()
                : Head(new TNode), Tail(Head.get()) {}

        TThreadsafeQueue(const TThreadsafeQueue& other) = delete;

        TThreadsafeQueue& operator=(const TThreadsafeQueue& other) = delete;

        bool TryPop(T& value) {
                return PopHead(value);
        }

        void Push(T&& newValue) {
                PushImpl(std::move(newValue));
        }

        void Push(const T& newValue) {
                PushImpl(newValue);
        }
};

class TFunctionWrapper {
private:
        struct TBaseImpl {
                virtual void Call() = 0;

                virtual ~TBaseImpl() {}
        };

        template<typename TFunctor>
        struct TImplType : TBaseImpl {
                TFunctor F;

                TImplType(TFunctor&& f)
                        : F(std::move(f)) {}

                void Call() override {
                        F();
                }
        };

        std::unique_ptr <TBaseImpl> ImplPtr;

public:

        template<typename TFunc>
        TFunctionWrapper(TFunc&& f)
                :   ImplPtr(new TImplType<TFunc>(std::forward<TFunc>(f))) {}

        void operator()() {
                ImplPtr->Call();
        }

        TFunctionWrapper() {}

        ~TFunctionWrapper() {}

        TFunctionWrapper(TFunctionWrapper&& functionWrapper)
                : ImplPtr(std::move(functionWrapper.ImplPtr)) {}

        TFunctionWrapper& operator=(TFunctionWrapper&& functionWrapper) {
                ImplPtr = std::move(functionWrapper.ImplPtr);
                return *this;
        }

        bool Empty() const {
                return ImplPtr == nullptr;
        }

private:
        TFunctionWrapper(const TFunctionWrapper&) = delete;

        TFunctionWrapper(TFunctionWrapper&) = delete;

        TFunctionWrapper& operator=(const TFunctionWrapper&) = delete;

        TFunctionWrapper& operator=(TFunctionWrapper&) = delete;
};

class TThreadJoiner {
private:
        std::vector <std::thread>& Threads;

public:
        TThreadJoiner(std::vector <std::thread>& threads)
                : Threads(threads) {}

        void Clear() {
                for (auto& thread : Threads) {
                        if (thread.joinable())
                                thread.join();
                }
        }


        ~TThreadJoiner() {
                Clear();
        }
};

class TThreadPool {
private:
        TThreadsafeQueue<TFunctionWrapper> TaskQueue;
        std::atomic_bool Done;
        std::vector <std::thread> Threads;
        TThreadJoiner ThreadJoiner;

        void Worker(
#ifdef __gnu_linux__
                uint32_t id
#endif
        ) {
#ifdef __gnu_linux__
                numa_set_interleave_mask(numa_all_nodes_ptr);
                PinToCore(id);
#endif
                while (!Done) {
                        TFunctionWrapper task;
                        if (TaskQueue.TryPop(task))
                                task();
//                        else
//                                std::this_thread::yield();
                }
        }

public:
        explicit TThreadPool(size_t threadsCount = 0)
                : ThreadJoiner(Threads) {
                Done = false;

                try {
                        Threads.reserve(threadsCount);
                        for (size_t i = 0; i < threadsCount; ++i)
                                Threads.push_back(std::thread(&TThreadPool::Worker, this
#ifdef __gnu_linux__
                                        , i + 1
#endif
                                ));
                }
                catch (...) {
                        Done = true;
                        throw;
                }
        }

        void Resize(size_t threadsCount) {
                Done = true;
                ThreadJoiner.Clear();
                Threads.clear();

                Done = false;
                Threads.reserve(threadsCount);
                for (size_t i = 0; i < threadsCount; ++i)
                        Threads.push_back(std::thread(&TThreadPool::Worker, this
#ifdef __gnu_linux__
                                , i + 1
#endif
                        ));

        }

        size_t NumThreads() const {
                return Threads.size();
        }

        ~TThreadPool() {
                Done = true;
        }

        template<typename TFunctor, typename... TArgs>
        std::future<typename std::result_of<TFunctor(TArgs...)>::type> Submit(TFunctor&& f, TArgs... args) {
                typedef typename std::result_of<TFunctor(TArgs...)>::type TResultType;
                std::packaged_task < TResultType() >
                task(std::bind(std::forward<TFunctor>(f), std::forward<TArgs>(args)...));
                std::future <TResultType> res(task.get_future());

                TaskQueue.Push(std::move(task));

                return res;
        }

        template<typename TFunctor, typename... TArgs>
        std::vector <std::future<typename std::result_of<TFunctor()>::type>> SubmitAll(TFunctor&& f, TArgs... args) {
                typedef typename std::result_of<TFunctor()>::type TResultType;

                std::vector <std::future<TResultType>> futures;
                futures.reserve(NumThreads());
                for (size_t i = 0; i < NumThreads(); ++i) {
                        std::packaged_task < TResultType() > task(std::bind(std::forward<TFunctor>(f),
                                                                            std::forward<TArgs>(args)...));

                        futures.push_back(task.get_future());
                        TaskQueue.Push(std::move(task));
                }

                return futures;
        }
};

class TThreadPoolWithTaskQueuePerThread {
private:
        using TQueue = CacheAlignedData<TThreadsafeQueue<TFunctionWrapper>>;

        std::atomic_bool Done;
        std::vector <std::thread> Threads;
        TThreadJoiner ThreadJoiner;
        std::unique_ptr<TQueue[]> Queues;

        void Worker(
#ifdef __gnu_linux__
                uint32_t core_id
#endif
        ) {
#ifdef __gnu_linux__
                PinToCore(core_id);
#endif
                while (!Done) {
                        TFunctionWrapper task;
                        if (Queues[core_id - 1].get().TryPop(task))
                                task();
//                        else
//                                std::this_thread::yield();
                }
#ifdef __gnu_linux__
                Unpin();
#endif
        }

public:
        explicit TThreadPoolWithTaskQueuePerThread(size_t threadsCount = 0)
                :       ThreadJoiner(Threads)
                ,       Queues(std::make_unique<TQueue[]>(threadsCount))
        {

                Done = false;

                try {
                        Threads.reserve(threadsCount);
                        for (size_t i = 0; i < threadsCount; ++i)
                                Threads.push_back(std::thread(&TThreadPoolWithTaskQueuePerThread::Worker, this
#ifdef __gnu_linux__
                                        , i + 1
#endif
                                ));
                }
                catch (...) {
                        Done = true;
                        throw;
                }
        }

        void Resize(size_t threadsCount) {
                Done = true;
                ThreadJoiner.Clear();
                Threads.clear();
                Queues = std::make_unique<TQueue[]>(threadsCount);

                Done = false;
                Threads.reserve(threadsCount);
                for (size_t i = 0; i < threadsCount; ++i)
                        Threads.push_back(std::thread(&TThreadPoolWithTaskQueuePerThread::Worker, this
#ifdef __gnu_linux__
                                , i + 1
#endif
                        ));

        }

        size_t NumThreads() const {
                return Threads.size();
        }

        void Clear() {
                Done = true;
                ThreadJoiner.Clear();
                Threads.clear();
        }

        ~TThreadPoolWithTaskQueuePerThread() {
                Done = true;
        }

        template<typename TFunctor, typename... TArgs>
        std::future<typename std::result_of<TFunctor(TArgs...)>::type> Submit(size_t thread_id, TFunctor&& f,
                                                                              TArgs... args) {
                typedef typename std::result_of<TFunctor(TArgs...)>::type TResultType;
                std::packaged_task < TResultType() >
                        task(std::bind(std::forward<TFunctor>(f), std::forward<TArgs>(args)...));
                std::future <TResultType> res(task.get_future());

                Queues[thread_id].get().Push(std::move(task));

                return res;
        }

        template<typename TFunctor, typename... TArgs>
        std::vector <std::future<typename std::result_of<TFunctor()>::type>> SubmitAll(size_t thread_id, TFunctor&& f,
                                                                                       TArgs... args) {
                typedef typename std::result_of<TFunctor()>::type TResultType;

                std::vector <std::future<TResultType>> futures;
                futures.reserve(NumThreads());
                for (size_t i = 0; i < NumThreads(); ++i) {
                        std::packaged_task < TResultType() > task(std::bind(std::forward<TFunctor>(f),
                                                                            std::forward<TArgs>(args)...));

                        futures.push_back(task.get_future());
                        Queues[thread_id].get().Push(std::move(task));
                }

                return futures;
        }
};

extern TThreadPoolWithTaskQueuePerThread g_thread_pool;

template<typename TFunctor>
static void submit_for_all(TFunctor functor) {
        std::vector<std::future<void>> futures;
        futures.reserve(g_thread_pool.NumThreads());
        for (uint32_t i = 0; i < g_thread_pool.NumThreads(); ++i) {
                if constexpr (function_traits<TFunctor>::arity == 0) {
                        futures.push_back(g_thread_pool.Submit(i, functor));
                }
                if constexpr (function_traits<TFunctor>::arity == 1) {
                        futures.push_back(g_thread_pool.Submit(i, functor, i + 1));
                }
        }
        if constexpr (function_traits<TFunctor>::arity == 0) {
                functor();
        }
        if constexpr (function_traits<TFunctor>::arity == 1) {
                functor(uint32_t(0));
        }

        std::for_each(futures.begin(), futures.end(), [&](auto& future) {
                future.get();
        });
};

template<typename TFunctor, typename TFunctorResult, typename TArg>
static typename std::result_of<TFunctorResult(TArg, TArg)>::type submit_for_all(TFunctor functor,
                                                                               TFunctorResult functor_result,
                                                                               const TArg& init_value) {
        std::vector<std::future<TArg>> futures;
        futures.reserve(g_thread_pool.NumThreads());
        for (uint32_t i = 0; i < g_thread_pool.NumThreads(); ++i) {
                if constexpr (function_traits<TFunctor>::arity == 0) {
                        futures.push_back(g_thread_pool.Submit(i, functor));
                }
                if constexpr (function_traits<TFunctor>::arity == 1) {
                        futures.push_back(g_thread_pool.Submit(i, functor, i + 1));
                }
        }

        using  res_type = typename std::result_of<TFunctorResult(TArg, TArg)>::type;

        res_type res = res_type();
        if constexpr (function_traits<TFunctor>::arity == 0) {
                res = functor_result(init_value, functor());
        }

        if constexpr (function_traits<TFunctor>::arity == 1) {
                res = functor_result(init_value, functor(uint32_t(0)));
        }

        std::for_each(futures.begin(), futures.end(), [&](auto& future) {
                res = functor_result(res, future.get());
        });
        return res;
};

template<typename TFunctor, typename TFunctorResult, typename TArg>
static void submit_for_all(TFunctor functor, TFunctorResult functor_result, TArg& result) {
        std::vector<std::future<TArg>> futures;
        futures.reserve(g_thread_pool.NumThreads());
        for (uint32_t i = 0; i < g_thread_pool.NumThreads(); ++i) {
                if constexpr (function_traits<TFunctor>::arity == 0) {
                        futures.push_back(g_thread_pool.Submit(i, functor));
                }
                if constexpr (function_traits<TFunctor>::arity == 1) {
                        futures.push_back(g_thread_pool.Submit(i, functor, i + 1));
                }
        }

        if constexpr (function_traits<TFunctor>::arity == 0) {
                functor_result(result, functor());
        }
        if constexpr (function_traits<TFunctor>::arity == 1) {
                functor_result(result, functor(uint32_t(0)));
        }

        std::for_each(futures.begin(), futures.end(), [&](auto& future) {
                functor_result(result, future.get());
        });
};

}