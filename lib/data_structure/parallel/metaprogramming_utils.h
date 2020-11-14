#pragma once

namespace parallel {
template <typename T>
struct function_traits : public function_traits<decltype(&T::operator())> {
};

template <typename class_type, typename _return_type, typename... args_types>
struct function_traits<_return_type(class_type::*)(args_types...) const> {
        static constexpr size_t arity = sizeof ... (args_types);
        using result_type = _return_type;
};

template <typename _return_type, typename... args_types>
struct function_traits<_return_type (*)(args_types...)> {
        static constexpr size_t arity = sizeof ... (args_types);
        using result_type = _return_type;
};
}
